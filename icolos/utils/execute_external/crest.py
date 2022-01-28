from icolos.utils.execute_external.execute import ExecutorBase
from icolos.utils.enums.program_parameters import CrestEnum


EE = CrestEnum()


class CrestExecutor(ExecutorBase):
    """For the execution of the "crest" binary."""

    def __init__(self, prefix_execution=None, binary_location=None):
        super().__init__(
            prefix_execution=prefix_execution, binary_location=binary_location
        )

    def execute(
        self, command: str, arguments: list, check=True, location=None, pipe_input=None
    ):
        # check, whether a proper executable is provided
        if command not in [EE.CREST]:
            raise ValueError(
                "Parameter command must be an dictionary of the internal Crest executable list."
            )

        return super().execute(
            command=command,
            arguments=arguments,
            check=check,
            location=location,
            pipe_input=pipe_input,
        )

    def is_available(self):
        try:
            result = self.execute(command=EE.CREST, arguments=[EE.CREST_H], check=True)

            if EE.CREST_HELP_IDENTIFICATION_STRING in result.stdout:
                return True
            return False
        except Exception as e:
            return False
