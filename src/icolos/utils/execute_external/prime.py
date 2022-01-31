from icolos.utils.enums.program_parameters import PrimeEnum, SchrodingerExecutablesEnum
from icolos.utils.execute_external.execute import ExecutorBase

SEE = SchrodingerExecutablesEnum()
EE = PrimeEnum()


class PrimeExecutor(ExecutorBase):
    """For the execution of the "prime_mmgbsa" binary."""

    def __init__(self, prefix_execution=None, binary_location=None):
        super().__init__(
            prefix_execution=prefix_execution, binary_location=binary_location
        )

    def execute(
        self, command: str, arguments: list, check=True, location=None, pipe_input=None
    ):
        # check, whether a proper executable is provided
        if command not in [EE.PRIME_MMGBSA]:
            raise ValueError(
                "Parameter command must be an dictionary of the internal Prime executable list."
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
            result = self.execute(
                command=EE.PRIME_MMGBSA, arguments=[EE.PRIME_HELP], check=True
            )

            if EE.PRIME_HELP_IDENTIFICATION_STRING in result.stdout:
                return True
            return False
        except Exception as e:
            return False
