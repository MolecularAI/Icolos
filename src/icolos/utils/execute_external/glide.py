from icolos.utils.enums.program_parameters import SchrodingerExecutablesEnum, GlideEnum
from icolos.utils.execute_external.execute import ExecutorBase

SEE = SchrodingerExecutablesEnum()
EE = GlideEnum()


class GlideExecutor(ExecutorBase):
    """For the execution of the "glide" binary."""

    def __init__(self, prefix_execution=None, binary_location=None):
        super().__init__(
            prefix_execution=prefix_execution, binary_location=binary_location
        )

    def execute(
        self, command: str, arguments: list, check=True, location=None, pipe_input=None
    ):
        # check, whether a proper executable is provided
        if command not in [EE.GLIDE]:
            raise ValueError(
                "Parameter command must be an dictionary of the internal Glide executable list."
            )

        # Note: It seems in former times, the call "glide" had to be changed to "$SCHRODINGER/glide" here.
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
                command=EE.GLIDE, arguments=[EE.GLIDE_HELP], check=True
            )

            if EE.GLIDE_HELP_IDENTIFICATION_STRING in result.stdout:
                return True
            return False
        except Exception as e:
            return False
