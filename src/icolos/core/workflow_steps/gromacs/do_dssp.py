from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum
from icolos.utils.enums.program_parameters import GromacsEnum
from icolos.utils.execute_external.gromacs import GromacsExecutor
from pydantic import BaseModel
from icolos.core.workflow_steps.step import _LE
from icolos.core.workflow_steps.gromacs.base import StepGromacsBase

_GE = GromacsEnum()
_SGE = StepGromacsEnum()
_SBE = StepBaseEnum


class StepGMXDoDSSP(StepGromacsBase, BaseModel):
    """
    Wrapper for gmx do_dssp binary, which in turn wraps dssp
    returns secondary structure data for the provided gromacs trajectory
    """

    def __init__(self, **data):
        super().__init__(**data)

        self._initialize_backend(executor=GromacsExecutor)
        self._check_backend_availability()

    def execute(self):
        tmp_dir = self._make_tmpdir()
        self.write_input_files(tmp_dir)

        structure_file = self.data.generic.get_argument_by_extension(_SGE.FIELD_KEY_TPR)
        traj_file = self.data.generic.get_argument_by_extension(_SGE.FIELD_KEY_XTC)
        arguments = self._parse_arguments(
            flag_dict={"-f": traj_file, "-s": structure_file, "-ssdump": "info.dat"}
        )

        pipe_input = (
            self.construct_pipe_arguments(
                tmp_dir, self.settings.additional[_SBE.PIPE_INPUT]
            )
            if _SBE.PIPE_INPUT in self.settings.additional.keys()
            and self.settings.additional[_SBE.PIPE_INPUT] is not None
            else None
        )

        result = self._backend_executor.execute(
            command=_GE.DO_DSSP,
            arguments=arguments,
            location=tmp_dir,
            pipe_input=pipe_input,
        )
        for line in result.stdout.split("\n"):
            self._logger_blank.log(line, _LE.DEBUG)
        self._logger.log(
            f"Completed execution for {self.step_id} successfully", _LE.INFO
        )

        self._parse_output(tmp_dir)
        self._remove_temporary(tmp_dir)
