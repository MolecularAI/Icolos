import os
import shutil
import tempfile

from icolos.utils.enums.program_parameters import TurbomoleEnum
from icolos.utils.execute_external.execute import ExecutorBase

EE = TurbomoleEnum()


class TurbomoleExecutor(ExecutorBase):
    """For the execution of the "turbomole" binary."""

    def __init__(self, prefix_execution=None, binary_location=None):
        super().__init__(
            prefix_execution=prefix_execution, binary_location=binary_location
        )

    def execute(
        self, command: str, arguments: list, check=True, location=None, pipe_input=None
    ):
        # check, whether a proper executable is provided
        if command not in [
            EE.TM_COSMOPREP,
            EE.TM_DEFINE,
            EE.TM_RIDFT,
            EE.TM_X2T,
            EE.TM_T2X,
            EE.CT_COSMOTHERM,
            EE.TM_JOBEX,
        ]:
            raise ValueError(
                "Parameter command must be an dictionary of the internal Turbomole executable list."
            )

        # TM accesses a folder specified in $TURBOTMPDIR to deposit the particular temporary files for a run; this is
        # system-wide, so parallel runs will interfere; also it is not removed automatically
        # TODO: find a more elegant solution; is this really necessary for all binaries or only "ridft" and "jobex"?
        tmp_dir = tempfile.mkdtemp()
        command = "".join(["export TURBOTMPDIR=", tmp_dir, " && ", command])

        result = super().execute(
            command=command,
            arguments=arguments,
            check=check,
            location=location,
            pipe_input=pipe_input,
        )

        if tmp_dir is not None and os.path.isdir(tmp_dir):
            shutil.rmtree(tmp_dir)
        return result

    def is_available(self):
        try:
            result = self.execute(command=EE.TM_RIDFT, arguments=[], check=True)

            if EE.TM_RIDFT_FAIL_IDENTIFICATION_STRING in result.stderr:
                return True
            return False
        except Exception as e:
            return False
