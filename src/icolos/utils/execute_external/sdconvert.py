from icolos.utils.enums.program_parameters import SchrodingerExecutablesEnum
from icolos.utils.execute_external.execute import ExecutorBase

SEE = SchrodingerExecutablesEnum()


class SDConvertExecutor(ExecutorBase):
    """For the execution of the "sdconvert" binary."""

    def __init__(self, prefix_execution=None, binary_location=None):
        super().__init__(
            prefix_execution=prefix_execution, binary_location=binary_location
        )

    def execute(
        self, command: str, arguments: list, check=True, location=None, pipe_input=None
    ):
        # check, whether a proper executable is provided
        if command not in [SEE.SDCONVERT]:
            raise ValueError(
                "Parameter command must be an dictionary of the internal sdconvert executable list."
            )

        # take care of the special path to "sdconvert"
        if command == SEE.SDCONVERT:
            command = SEE.SDCONVERT_CALL

        return super().execute(
            command=command,
            arguments=arguments,
            check=check,
            location=None,
            pipe_input=pipe_input,
        )

    def is_available(self):
        try:
            # raise NotImplementedError
            result = self.execute(
                command=SEE.SDCONVERT, arguments=SEE.SDCONVERT_HELP, check=False
            )

            if SEE.SDCONVERT_HELP_IDENTIFICATION_STRING in result.stdout:
                return True
            return False
        except Exception as e:
            return False
