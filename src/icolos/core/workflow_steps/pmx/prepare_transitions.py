from typing import Dict, List
from icolos.core.containers.perturbation_map import Edge
from icolos.core.workflow_steps.pmx.base import StepPMXBase
from icolos.core.workflow_steps.step import _LE
from pydantic import BaseModel
from icolos.utils.enums.program_parameters import (
    GromacsEnum,
    StepPMXEnum,
)
from icolos.utils.execute_external.pmx import PMXExecutor
from icolos.utils.general.parallelization import SubtaskContainer
import os


_PSE = StepPMXEnum()
_GE = GromacsEnum()


class StepPMXPrepareTransitions(StepPMXBase, BaseModel):
    """
    Prepare transitions: extract snapshots from equilibrium simulations, prepare .tpr files for each
    """

    def __init__(self, **data):
        super().__init__(**data)

        self._initialize_backend(executor=PMXExecutor)

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
        self._gromacs_executor.execute(
            _GE.TRJCONV, arguments=trjconv_args, pipe_input="echo System"
        )

        # move frame0.gro to frame80.gro
        cmd = "mv {0}/frame0.gro {0}/frame80.gro".format(tipath)
        os.system(cmd)

        self._clean_backup_files(tipath)

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
        self._extract_snapshots(eqpath, tipath)
        result = self._prepare_single_tpr(
            simpath=tipath,
            toppath=toppath,
            state=state,
            sim_type="transitions",
            framestart=1,
            framestop=81,
        )
        if result.returncode != 0:
            self._logger.log(f"WARNING, grompp has failed in {tipath}", _LE.WARNING)
            for line in result.stderr.split("\n"):
                self._logger.log(line, _LE.DEBUG)
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
        output_files = [
            f"ligand/stateA/run1/transitions/ti80.tpr",
            f"ligand/stateB/run1/transitions/ti80.tpr",
            f"complex/stateA/run1/transitions/ti80.tpr",
            f"complex/stateB/run1/transitions/ti80.tpr",
        ]
        results = []
        for subjob in batch:
            subjob_results = []
            for job in subjob:
                subjob_results.append(
                    all(
                        [
                            os.path.isfile(os.path.join(self.work_dir, job, f))
                            for f in output_files
                        ]
                    )
                )
            results.append(subjob_results)
        return results
