from typing import List
from icolos.core.workflow_steps.pmx.base import StepPMXBase
from icolos.core.workflow_steps.step import _LE
from pydantic import BaseModel
from icolos.utils.enums.program_parameters import (
    GromacsEnum,
)
from icolos.utils.execute_external.gromacs import GromacsExecutor
from icolos.utils.general.parallelization import SubtaskContainer
import os

_GE = GromacsEnum()


class StepPMXPrepareTransitions(StepPMXBase, BaseModel):
    """
    Prepare transitions: extract snapshots from equilibrium simulations, prepare .tpr files for each
    """

    def __init__(self, **data):
        super().__init__(**data)

        self._initialize_backend(executor=GromacsExecutor)

    def execute(self):
        edges = [e.get_edge_id() for e in self.get_edges()]

        self.execution.parallelization.max_length_sublists = 1
        self._subtask_container = SubtaskContainer(
            max_tries=self.execution.failure_policy.n_tries
        )
        self._subtask_container.load_data(edges)
        self._execute_pmx_step_parallel(
            run_func=self.prepare_transitions,
            step_id="pmx prepare_transitions",
            result_checker=self._check_result,
        )

    def _extract_snapshots(self, eqpath, tipath):
        tpr = "{0}/tpr.tpr".format(eqpath)
        trr = "{0}/traj.trr".format(eqpath)
        frame = "{0}/frame.gro".format(tipath)

        trjconv_args = {
            "-s": tpr,
            "-f": trr,
            "-o": frame,
            "-sep": "",
            "-ur": "compact",
            "-pbc": "mol",
            "-b": 2000,
        }
        trjconv_args = self.get_arguments(trjconv_args)
        self._backend_executor.execute(
            _GE.TRJCONV, arguments=trjconv_args, pipe_input="echo System", check=False
        )

        last_frame = len([f for f in os.listdir(tipath) if f.startswith("frame")])
        self._logger.log(f"Extracted {last_frame} frames", _LE.DEBUG)

        self._clean_backup_files(tipath)
        # once frames are extracted, remove the large trr file from the equilibrium run
        files_to_remove = [trr, tpr, os.path.join(eqpath, "frame.gro")]
        for f in files_to_remove:
            os.remove(f)

    def _prepare_system(self, edge: str, state: str, wp: str, r: int, toppath: str):
        eqpath = self._get_specific_path(
            workPath=self.work_dir,
            edge=edge,
            wp=wp,
            state=state,
            r=r,
            sim="eq",
        )
        tipath = self._get_specific_path(
            workPath=self.work_dir,
            edge=edge,
            wp=wp,
            state=state,
            r=r,
            sim="transitions",
        )
        # if the trr file exists and snapshots have not been extracted in a previous run
        if os.path.isfile(os.path.join(eqpath, "traj.trr")) and not os.path.isfile(
            os.path.join(tipath, "frame0.gro")
        ):
            self._extract_snapshots(eqpath, tipath)
        else:
            self._logger.log("Skipping frame extraction, already present")

        self._prepare_single_tpr(
            simpath=tipath,
            toppath=toppath,
            state=state,
            sim_type="transitions",
            executor=self._backend_executor,
        )

        self._clean_backup_files(tipath)

    def prepare_transitions(self, jobs: List[str]):
        for edge in jobs:
            ligTopPath = self._get_specific_path(
                workPath=self.work_dir, edge=edge, wp="ligand"
            )
            protTopPath = self._get_specific_path(
                workPath=self.work_dir, edge=edge, wp="complex"
            )
            for state in self.states:
                for r in range(1, self.get_perturbation_map().replicas + 1):

                    self._logger.log(
                        f"Preparing transitions: {edge}, {state}, run {r}", _LE.DEBUG
                    )
                    self._prepare_system(
                        edge=edge, state=state, wp="ligand", r=r, toppath=ligTopPath
                    )
                    self._prepare_system(
                        edge=edge, state=state, wp="complex", r=r, toppath=protTopPath
                    )

    def _check_result(self, batch: List[List[str]]) -> List[List[bool]]:
        """
        Look in each hybridStrTop dir and check the output pdb files exist for the edges
        """
        replicas = self.get_perturbation_map().replicas
        output_paths = []
        for i in range(1, replicas + 1):
            output_paths.append(f"ligand/stateA/run{i}/transitions")
            output_paths.append(f"ligand/stateB/run{i}/transitions")
            output_paths.append(f"complex/stateA/run{i}/transitions")
            output_paths.append(f"complex/stateB/run{i}/transitions")

        results = []
        for subjob in batch:
            subjob_results = []
            for job in subjob:
                # check number of frames extracted == number tpr files for each leg

                num_frames = [
                    len(
                        [
                            f
                            for f in os.listdir(os.path.join(self.work_dir, job, f))
                            if f.startswith("frame")
                        ]
                    )
                    for f in output_paths
                ]
                num_tprs = [
                    len(
                        [
                            f
                            for f in os.listdir(os.path.join(self.work_dir, job, f))
                            if f.startswith("ti")
                        ]
                    )
                    for f in output_paths
                ]
                # confirms that frames have been extracted, and we have a tpr file generated for each ti frame
                subjob_results.append(
                    num_tprs == num_frames and all(n > 0 for n in num_frames)
                )
            results.append(subjob_results)
        return results
