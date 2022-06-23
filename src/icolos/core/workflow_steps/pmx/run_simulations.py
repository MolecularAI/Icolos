from typing import List
from icolos.core.workflow_steps.pmx.base import StepPMXBase
from pydantic import BaseModel
from icolos.core.workflow_steps.step import _LE
import numpy as np
from icolos.utils.enums.execution_enums import ExecutionPlatformEnum
from icolos.utils.enums.program_parameters import (
    StepPMXEnum,
)
from icolos.utils.execute_external.slurm_executor import SlurmExecutor
from icolos.utils.general.parallelization import Subtask, SubtaskContainer
import os

_PSE = StepPMXEnum()
_EPE = ExecutionPlatformEnum


class StepPMXRunSimulations(StepPMXBase, BaseModel):
    """
    Runs md simulations, unwraps into pool of individual jobs, parallelized over available GPUs
    """

    sim_type: str = ""

    def __init__(self, **data):
        super().__init__(**data)

        # Note: if youre running the job on, for example, a workstation, without slurm, this will simply execute the scripts directly (the slurm header is simply ignored in this case)
        self._initialize_backend(executor=SlurmExecutor)

    def execute(self):

        if self.run_type == "rbfe":
            edges = [e.get_edge_id() for e in self.get_edges()]
        elif self.run_type == "abfe":
            edges = [c.get_index_string() for c in self.get_compounds()]
        self.sim_type = self._get_additional_setting(_PSE.SIM_TYPE)
        assert (
            self.sim_type in self.mdp_prefixes.keys()
        ), f"sim type {self.sim_type} not recognised!"

        # prepare and pool jobscripts, unroll replicas,  etc
        job_pool = self._prepare_job_pool(edges)
        self._logger.log(
            f"Prepared {len(job_pool)} jobs for {self.sim_type} simulations",
            _LE.DEBUG,
        )

        self._subtask_container = SubtaskContainer(
            max_tries=self.execution.failure_policy.n_tries
        )
        self._subtask_container.load_data(job_pool)
        result_checker = (
            self._inspect_dhdl_files
            if self.sim_type == "transitions"
            else self._inspect_log_files
        )
        # simulations are run without the normal parallelizer
        # this relies upon job IDs from Slurm
        if self.execution.platform == _EPE.SLURM:
            self._run_job_pool(run_func=self._run_single_job)
        else:
            # simulations are run locally
            self._execute_pmx_step_parallel(
                run_func=self._execute_command,
                step_id="pmx_run_simulations",
                # for run_simulations, because batch efficiency is crucial, we do this prior to batching
                prune_completed=False,
                result_checker=result_checker,
            )
        # clean up large files from completed transition jobs only
        # if self.sim_type == "transitions":
        #     self._logger.log("Cleaning up transition output files", _LE.DEBUG)
        #     to_clean_list = []
        #     for job in job_pool:
        #         # search branch, state, run
        #         to_clean_list.append(
        #             glob(f"{self.work_dir}/{job}/*/*/*/transitions/frame*.gro")
        #         )
        #     for file in to_clean_list:
        #         os.remove(file)

    def get_mdrun_command(
        self,
        tpr=None,
        ener=None,
        confout=None,
        mdlog=None,
        trr=None,
    ):

        mdrun_binary = self._get_additional_setting(
            _PSE.MDRUN_EXECUTABLE, default="gmx mdrun"
        )
        if self.sim_type in ("em", "eq", "npt", "nvt"):

            job_command = [
                mdrun_binary,
                "-s",
                tpr,
                "-e",
                ener,
                "-c",
                confout,
                "-o",
                trr,
                "-g",
                mdlog,
            ]
            for flag in self.settings.arguments.flags:
                job_command.append(str(flag))
            for key, value in self.settings.arguments.parameters.items():
                job_command.append(str(key))
                job_command.append(str(value))

        elif self.sim_type == "transitions":
            # need to add many job commands to the slurm file, one for each transition
            sim_path = os.path.dirname(tpr)
            tpr_files = [f for f in os.listdir(sim_path) if f.endswith("tpr")]
            job_command = []
            # note, these will not be returned in numerical order!
            for file in tpr_files:
                # grab the index in the tpr file
                tpr_idx = file.split(".tpr")[0][2:]
                dhdl_file = os.path.join(os.path.dirname(mdlog), f"dhdl{tpr_idx}.xvg")
                if not os.path.isfile(dhdl_file):
                    single_command = [
                        mdrun_binary,
                        "-s",
                        file,
                        "-e",
                        ener,
                        "-c",
                        confout,
                        "-dhdl",
                        f"dhdl{tpr_idx}.xvg",
                        "-o",
                        trr,
                        "-g",
                        mdlog,
                    ]
                    for flag in self.settings.arguments.flags:
                        single_command.append(str(flag))
                    for key, value in self.settings.arguments.parameters.items():
                        single_command.append(str(key))
                        single_command.append(str(value))

                    single_command.append("\n")
                    single_command.append(f"\nrm {os.path.dirname(ener)}/*#\n")
                    job_command += single_command
                else:
                    self._logger.log(
                        f"dhdl file for transition {tpr_idx} in {os.path.dirname(mdlog)} already exists, skipping",
                        _LE.DEBUG,
                    )

        return job_command

    def _prepare_single_job(self, edge: str, wp: str, state: str, r: int):
        """
        Construct a slurm job file in that job's directory, return the path to the batch script
        """
        simpath = self._get_specific_path(
            workPath=self.work_dir,
            edge=edge,
            wp=wp,
            state=state,
            r=r,
            sim=self.sim_type,
        )
        tpr = "{0}/tpr.tpr".format(simpath)
        ener = "{0}/ener.edr".format(simpath)
        confout = "{0}/confout.gro".format(simpath)
        mdlog = "{0}/md.log".format(simpath)
        trr = "{0}/traj.trr".format(simpath)
        if self.sim_type != "transitions":
            try:
                with open(os.path.join(simpath, "md.log"), "r") as f:
                    lines = f.readlines()
                sim_complete = any(["Finished mdrun" in l for l in lines])
            except FileNotFoundError:
                sim_complete = False

        # handle transitions
        else:
            # cannot reliably check that all sims for all edges have completed here, this will be checked in get_mdrun_command which will skip completed perturbations if dhdl exists
            sim_complete = False

        if not sim_complete:
            self._logger.log(
                f"Preparing: {wp} {edge} {state} run{r}, simType {self.sim_type}",
                _LE.DEBUG,
            )
            job_command = self.get_mdrun_command(
                tpr=tpr,
                trr=trr,
                ener=ener,
                confout=confout,
                mdlog=mdlog,
            )
            # empty list indicates all transitions completed
            if job_command:
                job_command = " ".join(job_command)
                batch_file = self._backend_executor.prepare_batch_script(
                    job_command, arguments=[], location=simpath
                )
                return os.path.join(simpath, batch_file)

        return

    def _prepare_job_pool(self, edges: List[str]):
        replicas = (
            self.get_perturbation_map().replicas
            if self.get_perturbation_map() is not None
            else 1
        )
        batch_script_paths = []
        # load in the protein jobs first, queue is FIFO
        for branch in self.therm_cycle_branches:
            for edge in edges:
                for r in range(1, replicas + 1):
                    for state in self.states:
                        path = self._prepare_single_job(
                            edge=edge, wp=branch, state=state, r=r
                        )
                        if path is not None:
                            batch_script_paths.append(path)
        return batch_script_paths

    def _run_single_job(self, job: str):
        """
        Execute the simulation for a single batch script
        """
        self._logger.log(
            f"Starting execution of job {job.split('/')[-1]}.",
            _LE.DEBUG,
        )
        location = os.path.dirname(job)
        job_id = self._backend_executor.execute(
            tmpfile=job, location=location, check=False, wait=False
        )
        return job_id

    def _execute_command(self, jobs: List[Subtask]):
        for idx, job in enumerate(jobs):
            self._logger.log(
                f"Starting execution of job {job.split('/')[-1]}. ",
                _LE.DEBUG,
            )
            self._logger.log(f"Batch progress: {idx+1}/{len(jobs)}", _LE.DEBUG)
            location = os.path.dirname(job)
            result = self._backend_executor.execute(
                tmpfile=job, location=location, check=False
            )
            self._logger.log(
                f"Execution for job {job} completed with status: {result}", _LE.DEBUG
            )

    def _inspect_log_files(self, jobs: List[str]) -> List[List[bool]]:
        """
        Check the md.log files in the edge's job dir, return
        :param jobs: list of paths to the batch scripts
        """
        results = []
        for subtask in jobs:
            subtask_results = []
            for sim in subtask:
                location = os.path.join(os.path.dirname(sim), "md.log")
                if os.path.isfile(location):
                    with open(location, "r") as f:
                        lines = f.readlines()
                    subtask_results.append(
                        any(["Finished mdrun" in l for l in lines[-20:]])
                    )
                else:
                    subtask_results.append(False)
            results.append(subtask_results)
        return results

    def _inspect_dhdl_files(self, jobs: List[str]) -> List[List[bool]]:
        # check there is a dhdl file for each tpr file
        results = []
        for subtask in jobs:
            subtask_results = []
            for sim in subtask:

                location = os.path.join(os.path.dirname(sim))
                n_snapshots = len(
                    [f for f in os.listdir(location) if f.endswith("tpr")]
                )
                dhdl_files = [f"dhdl{i}.xvg" for i in range(1, n_snapshots + 1)]

                if not all(
                    [os.path.isfile(os.path.join(location, f)) for f in dhdl_files]
                ):
                    subtask_results.append(False)
                else:
                    subtask_results.append(True)

            results.append(subtask_results)
        return results
