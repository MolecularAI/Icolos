from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum
from icolos.utils.enums.program_parameters import GromacsEnum
from icolos.core.workflow_steps.gromacs.base import StepGromacsBase
from icolos.utils.execute_external.gromacs import GromacsExecutor
from pydantic import BaseModel
from icolos.core.workflow_steps.step import _LE

_GE = GromacsEnum()
_SGE = StepGromacsEnum()
_SBE = StepBaseEnum


class StepGMXTrjconv(StepGromacsBase, BaseModel):
    """
    Postprocessing step for gromacs trajectories
    Mostly used for removing pbc, fitting trajectory etc.
    """

    def __init__(self, **data):
        super().__init__(**data)

        self._initialize_backend(executor=GromacsExecutor)
        self._check_backend_availability()

    def execute(self):

        tmp_dir = self._make_tmpdir()
        self._write_input_files(tmp_dir)

        xtc_file = self.data.generic.get_argument_by_extension(_SGE.FIELD_KEY_XTC)
        flag_dict = {
            "-s": self.data.generic.get_argument_by_extension(_SGE.FIELD_KEY_TPR),
            "-f": xtc_file,
            "-o": xtc_file,
        }

        arguments = self._parse_arguments(flag_dict=flag_dict)
        result = self._backend_executor.execute(
            command=_GE.TRJCONV,
            arguments=arguments,
            location=tmp_dir,
            check=True,
            pipe_input=self.construct_pipe_arguments(
                tmp_dir, self.settings.additional[_SBE.PIPE_INPUT]
            ),
        )
        for line in result.stdout.split("\n"):
            self._logger_blank.log(line, _LE.DEBUG)
        self._parse_output(tmp_dir)
        self._remove_temporary(tmp_dir)
