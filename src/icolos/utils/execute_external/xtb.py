from icolos.utils.enums.program_parameters import XTBEnum
from icolos.utils.execute_external.execute import ExecutorBase

EE = XTBEnum()


class XTBExecutor(ExecutorBase):
    """For the execution of the "xtb" binary."""

    def __init__(self, prefix_execution=None, binary_location=None):
        super().__init__(
            prefix_execution=prefix_execution, binary_location=binary_location
        )

    def execute(
        self, command: str, arguments: list, check=True, location=None, pipe_input=None
    ):
        # check, whether a proper executable is provided
        if command not in [EE.XTB]:
            raise ValueError(
                "Parameter command must be an dictionary of the internal XTB executable list."
            )

        return super().execute(
            command=command,
            arguments=arguments,
            check=check,
            location=None,
            pipe_input=pipe_input,
        )

    def is_available(self):
        try:
            result = self.execute(command=EE.XTB, arguments=[EE.XTB_HELP], check=True)

            if EE.XTB_HELP_IDENTIFICATION_STRING in result.stderr:
                return True
            return False
        except Exception as e:
            return False
