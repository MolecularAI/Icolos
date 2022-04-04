import tempfile
from typing import List
from icolos.core.containers.gmx_state import GromacsState
from icolos.utils.enums.execution_enums import ExecutionPlatformEnum
from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum
from icolos.utils.enums.program_parameters import GromacsEnum
from icolos.core.workflow_steps.gromacs.base import StepGromacsBase
from pydantic import BaseModel
from icolos.core.workflow_steps.step import _LE
import os

from icolos.utils.execute_external.gromacs import GromacsExecutor
from icolos.utils.general.parallelization import Parallelizer, SubtaskContainer

_GE = GromacsEnum()
_SGE = StepGromacsEnum()
_SBE = StepBaseEnum
_ERE = ExecutionPlatformEnum


class StepGMXMDrun(StepGromacsBase, BaseModel):
    """
    Launch gmx mdrun
    """

    topol: GromacsState = None

    def __init__(self, **data):
        super().__init__(**data)

        self._initialize_backend(executor=GromacsExecutor)
        self._check_backend_availability()

    def _get_log_file(self, tmp_dir):
        """
        Find and parse the log file
        """
        log_file = [f for f in os.listdir(tmp_dir) if f.endswith(".log")]
        assert len(log_file) == 1
        with open(os.path.join(tmp_dir, log_file[0]), "r") as f:
            data = f.readlines()
        return data

    def _tail_log_file(self, tmp_dir):
        """
        Log the last 50 lines of the log file to capture performance metrics from the run

        """
        log_file = self._get_log_file(tmp_dir)

        for line in log_file[-50:]:
            self._logger_blank.log(line, _LE.INFO)

    def execute_mdrun(self, path: str, index: int):
        """
        Make a single call to mdrun
        """
        flag_dict = (
            {
                "-s": _SGE.STD_TPR,
                "-c": _SGE.STD_STRUCTURE,
                "-x": _SGE.STD_XTC,
            }
            if not self.data.generic.get_files_by_extension("cpt")
            else {"-cpi", self.data.generic.get_argument_by_extension("cpt")}
        )

        arguments = self._parse_arguments(flag_dict)
        self._backend_executor.execute(
            command=_GE.MDRUN, arguments=arguments, location=path, check=True
        )

    def execute_parallel_simulations(self, work_dirs):
        # attach the index of the workdir
        work_dirs = [(idx, wkdir) for idx, wkdir in enumerate(work_dirs)]
        self._subtask_container = SubtaskContainer(
            max_tries=self.execution.failure_policy.n_tries
        )
        self._subtask_container.load_data(work_dirs)
        parallelizer = Parallelizer(func=self.execute_mdrun)
        n = 1
        while self._subtask_container.done() is False:
            next_batch = self._get_sublists(get_first_n_lists=self._get_number_cores())
            _ = [sub.increment_tries() for element in next_batch for sub in element]
            _ = [sub.set_status_failed() for element in next_batch for sub in element]
            paths, indices = self.prepare_jobs(next_batch)
            parallelizer.execute_parallel(path=paths, index=indices)
            n += 1

    def prepare_jobs(self, batch) -> List[tuple]:
        paths, indices = [], []
        for task in batch:
            for element in task:
                # tuple of (idx, dirpath)
                paths.append(element.data[1])
                indices.append(element.data[0])
        return paths, indices

    def run_single_tpr(self, tmp_dir: str):
        """
        Normal gmx mdrun call, if multiple structures are loaded into the topology, run them in parallel according to the parallelizer settings
        """
        # if we have multiple structures, run the simulations externally, in parallel
        work_dirs = [tempfile.mkdtemp(dir=tmp_dir) for _ in range(len(self.topol.tprs))]

        # prepare tmpdirs with tpr files
        for path, tpr in zip(work_dirs, self.topol.tprs.values()):
            tpr.write(path)

        # if > 1, instantiate a parallelizer, load the paths in and execute in parallel, user should be using the slurm/SGE interface to request extern resources
        if len(work_dirs) > 1:
            self.execute_parallel_simulations(work_dirs)
        else:
            tmp_dir = work_dirs[0]
            self.execute_mdrun(tmp_dir, index=0)

        # now parse the outputs
        for index, path in enumerate(work_dirs):
            # set a structure other than confout.gro e.g. if a pdb output has been set
            struct = (
                self.settings.arguments.parameters["-c"]
                if "-c" in self.settings.arguments.parameters.keys()
                else _SGE.STD_STRUCTURE
            )
            self.topol.set_structure(path, file=struct, index=index)
            self.topol.set_trajectory(path, index=index)
            self.topol.set_log(path, index=index)

    def run_multidir_sim(self, tmp_dir: str):
        """
        Runs a multidir simulation, allowing for replex simulations.  Several conditions are required for this running mode
        1) the previous step in the workflow should have been an iterator to produce n tpr files.  This must have been run with single_dir mode ON and remove_temprorary_files OFF, so we can extract files from those workflows' tmpdirs
        """
        if not self.execution.platform == _ERE.SLURM:
            self._logger.log(
                "WARNING: Running HREX simulation without external resources! Normally this should be run as a separate batch job",
                _LE.WARNING,
            )

        # extract the tprs from the topol object, write to separate tmpdirs
        work_dirs = [tempfile.mkdtemp(dir=tmp_dir) for _ in range(len(self.topol.tprs))]
        self._logger.log(
            f"Initiating gmx multidir run in directories {', '.join(work_dirs)}",
            _LE.DEBUG,
        )
        for path, tpr in zip(work_dirs, self.topol.tprs.values()):
            tpr.write(path)

        # note, this must be a multiple of the number of simulations
        tasks = self.execution.resources.tasks
        # map the PP and PME tasks to the GPUs

        command = f"mpirun -np {tasks} gmx_mpi mdrun -multidir {' '.join(work_dirs)}"
        arguments = self._parse_arguments(flag_dict={"-x": _SGE.STD_XTC})
        self._backend_executor.execute(
            command=command, arguments=arguments, location=tmp_dir, check=True
        )
        # udpate the structures to the new coordinates
        for i, work_dir in enumerate(work_dirs):
            self.topol.set_structure(work_dir, index=i)
            self.topol.set_trajectory(work_dir, index=i)
            self.topol.set_tpr(work_dir, index=i)
            self.topol.set_log(work_dir, index=i)

    def execute(self):

        tmp_dir = self._make_tmpdir()
        self.topol = self.get_topol()
        self.execution.parallelization.max_length_sublists = 1
        # pickle the topol to the mdrun dir, if something goes wrong/the job dies, the workflow can be picked up where we left off by unpickling the topology object
        self.pickle_topol(self.topol, tmp_dir)
        multidir = self.get_additional_setting(_SGE.MULTIDIR, default=False)
        if multidir:
            self.run_multidir_sim(tmp_dir)
        else:
            self.run_single_tpr(tmp_dir)
        self._remove_temporary(tmp_dir)
