from icolos.utils.enums.program_parameters import GromacsEnum
from icolos.core.workflow_steps.gromacs.base import StepGromacsBase
from icolos.utils.execute_external.gromacs import GromacsExecutor
from pydantic import BaseModel
from icolos.core.workflow_steps.step import _LE
import os

_GE = GromacsEnum()


class StepGMXselect(StepGromacsBase, BaseModel):
    """
    Run gromacs rmsd calculation on trajectory
    """

    def __init__(self, **data):
        super().__init__(**data)

        self._initialize_backend(executor=GromacsExecutor)
        self._check_backend_availability()

    def execute(self):

        tmp_dir = self._make_tmpdir()

        # write out generic files
        self.write_input_files(tmp_dir)

        flag_dict = {
            "-s": self.data.generic.get_argument_by_extension("tpr"),
            "-f": self.data.generic.get_argument_by_extension("xtc"),
            "-on": "index.ndx",
            "-select": self.settings.additional["selection"],
        }

        arguments = self._parse_arguments(flag_dict=flag_dict)
        result = self._backend_executor.execute(
            command=_GE.SELECT,
            arguments=arguments,
            location=tmp_dir,
            check=True,
        )

        self._log_result(result)
        self._parse_output(tmp_dir)
        self._remove_temporary(tmp_dir)
