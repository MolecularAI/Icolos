import tempfile
from icolos.core.containers.gromacs_topol import GromacsTopol
from icolos.utils.enums.execution_enums import ExecutionPlatformEnum
from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum
from icolos.utils.enums.program_parameters import GromacsEnum
from icolos.core.workflow_steps.gromacs.base import StepGromacsBase
from pydantic import BaseModel
from icolos.core.workflow_steps.step import _LE
import os

from icolos.utils.execute_external.gromacs import GromacsExecutor

_GE = GromacsEnum()
_SGE = StepGromacsEnum()
_SBE = StepBaseEnum
_ERE = ExecutionPlatformEnum


class StepGMXMDrun(StepGromacsBase, BaseModel):
    """
    Launch gmx mdrun
    """

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

    def run_single_tpr(self, tmp_dir: str, topol: GromacsTopol):
        # if we're simulating a protein, we need to modify the topol file to include the correct index groups \
        # to allow ligand restraint.  This means an ndx file must be specified in the json
        self.write_input_files(tmp_dir)
        # append _out to the xtc file name
        xtc_output_file = self.generate_output_file(_SGE.STD_XTC)
        arguments = self._parse_arguments(
            flag_dict={
                "-s": self.data.generic.get_argument_by_extension(_SGE.FIELD_KEY_TPR),
                "-c": _SGE.STD_STRUCTURE,
                "-x": xtc_output_file,
            }
        )
        self._backend_executor.execute(
            command=_GE.MDRUN, arguments=arguments, location=tmp_dir, check=True
        )

        self._tail_log_file(tmp_dir)
        self._logger.log(
            f"Completed execution for {self.step_id} successfully", _LE.INFO
        )
        self._parse_output(tmp_dir)
        if "-c" in self.settings.arguments.parameters.keys():
            topol.set_structure(tmp_dir, self.settings.arguments.parameters["-c"])
        else:
            topol.set_structure(tmp_dir)

    def run_multidir_sim(self, tmp_dir: str, topol: GromacsTopol):
        """
        Runs a multidir simulation, allowing for replex simulations.  Several conditions are required for this running mode
        1) the previous step in the workflow should have been an iterator to produce n tpr files.  This must have been run with single_dir mode ON and remove_temprorary_files OFF, so we can extract files from those workflows' tmpdirs

        """
        if not self.execution.platform == _ERE.SLURM:
            self._logger.log(
                "WARNING: Running HREX simulation using workflow's resources! Normally this should be run as a separate batch job",
                _LE.WARNING,
            )

        # extract the tprs from the topol object, write to separate tmpdirs
        work_dirs = [
            tempfile.mkdtemp(dir=tmp_dir) for _ in range(len(topol.structures))
        ]
        self._logger.log(
            f"Initiating gmx multidir run in directories {', '.join(work_dirs)}",
            _LE.DEBUG,
        )
        for path, tpr in zip(work_dirs, topol.tprs):
            tpr.write(path)

        # note, this must be a multiple of the number of simulations
        threads = self.execution.resources.other_args["--threads"]
        # map the PP and PME tasks to the GPUs

        command = f"mpirun -np {threads} gmx_mpi mdrun -multidir {' '.join(work_dirs)}"
        arguments = self._parse_arguments(flag_dict={"-x": _SGE.STD_XTC})
        self._backend_executor.execute(
            command=command, arguments=arguments, location=tmp_dir, check=True
        )
        # udpate the structures to the new coordinates
        for i, work_dir in enumerate(work_dirs):
            topol.set_structure(work_dir, index=i)
            topol.set_trajectory(work_dir, index=i)
            topol.set_tpr(work_dir, index=i)

    def execute(self):

        tmp_dir = self._make_tmpdir()
        topol = self.get_topol()
        multidir = self.get_additional_setting(_SGE.MULTIDIR, default=False)
        if multidir:
            self.run_multidir_sim(tmp_dir, topol=topol)
        else:
            self.run_single_tpr(tmp_dir, topol=topol)
        self._remove_temporary(tmp_dir)
