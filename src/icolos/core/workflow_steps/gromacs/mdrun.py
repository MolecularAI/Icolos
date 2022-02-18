from icolos.core.containers.gromacs_topol import GromacsTopol
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

    def run_multidir_sim(self, tmp_dir, topol):
        """
        Runs a multidir simulation, allowing for replex simulations.  Several conditions are required for this running mode
        1) the previous step in the workflow should have been an iterator to produce n tpr files.  This must have been run with single_dir mode ON and remove_temprorary_files OFF, so we can extract files from those workflows' tmpdirs

        """
        dispatcher_step = [
            s
            for s in self.get_workflow_object()._initialized_steps
            if s.type == _SBE.DISPATCHER
        ][-1]
        workflows = dispatcher_step.workflows
        # now we have the workdirs which should contain the final tprs from the production grompp step
        work_dirs = [wf.work_dir for wf in workflows]
        # note, this must be a multiple of the number of simulations
        threads = self.execution.resources.other_args["threads"]
        n_gpus = int(self.execution.resources.gres.split(":")[-1])
        # map the PP and PME tasks to the GPUs
        n_sims = len(work_dirs)
        gputask_str = ""
        for _ in n_sims:
            gputask_str += "01"
        command = f"mpirun -np {threads} gmx_mpi mdrun -multidir {' '.join(work_dirs)} --gputasks 0000011111"

    def execute(self):

        tmp_dir = self._make_tmpdir()
        topol = self.get_topol()
        run_mode = self.get_additional_setting(_SGE.RUN_MODE, default=_SGE.SERIAL)
        if run_mode == _SGE.SERIAL:
            self.run_single_tpr(tmp_dir, topol=topol)
        elif run_mode == _SGE.MULTIDIR:
            self.run_multidir_sim(tmp_dir)
        self._remove_temporary(tmp_dir)
