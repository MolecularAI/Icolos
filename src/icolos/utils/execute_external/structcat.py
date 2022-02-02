from icolos.utils.enums.program_parameters import SchrodingerExecutablesEnum
from icolos.utils.execute_external.execute import ExecutorBase

SEE = SchrodingerExecutablesEnum()


class StructcatExecutor(ExecutorBase):
    """For the execution of the "structcat" binary."""

    def __init__(self, prefix_execution=None, binary_location=None):
        super().__init__(
            prefix_execution=prefix_execution, binary_location=binary_location
        )

    def execute(
        self, command: str, arguments: list, check=True, location=None, pipe_input=None
    ):
        # check, whether a proper executable is provided
        if command not in [SEE.STRUCTCAT]:
            raise ValueError(
                "Parameter command must be an dictionary of the internal structcat executable list."
            )

        # take care of the special path to "structcat"
        if command == SEE.STRUCTCAT:
            command = SEE.STRUCTCAT_CALL

        return super().execute(
            command=command,
            arguments=arguments,
            check=check,
            location=None,
            pipe_input=pipe_input,
        )

    def is_available(self):
        try:
            result = self.execute(
                command=SEE.STRUCTCAT, arguments=SEE.STRUCTCAT_HELP, check=False
            )

            if SEE.STRUCTCAT_HELP_IDENTIFICATION_STRING in result.stdout:
                return True
            return False
        except Exception as e:
            return False
