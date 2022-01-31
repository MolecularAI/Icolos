from icolos.utils.enums.step_enums import StepGromacsEnum
from icolos.utils.enums.program_parameters import GromacsEnum
from icolos.core.workflow_steps.gromacs.base import StepGromacsBase
from pydantic import BaseModel
from icolos.core.workflow_steps.step import _LE
import os

_GE = GromacsEnum()
_SGE = StepGromacsEnum()


class StepGMXMDrun(StepGromacsBase, BaseModel):
    """
    Launch gmx mdrun
    """

    def __init__(self, **data):
        super().__init__(**data)

        self._initialize_backend(executor=self._get_gromacs_executor())
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

    def execute(self):

        tmp_dir = self._make_tmpdir()
        # if we're simulating a protein, we need to modify the topol file to include the correct index groups \
        # to allow ligand restraint.  This means an ndx file must be specified in the json
        self._write_input_files(tmp_dir)
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
        self._remove_temporary(tmp_dir)
