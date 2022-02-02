from icolos.utils.execute_external.execute import ExecutorBase


class CressetExecutor(ExecutorBase):
    """For the execution of Cresset binaries binary."""

    def __init__(self, prefix_execution=None, binary_location=None):
        super().__init__(
            prefix_execution=prefix_execution, binary_location=binary_location
        )

    def execute(
        self, command: str, arguments: list, check=True, location=None, pipe_input=None
    ):
        # check, whether a proper executable is provided
        # if command not in [EE.OMEGA]:
        #     raise ValueError(
        #         "Parameter command must be an dictionary of the internal Omega executable list."
        #     )

        return super().execute(
            command=command,
            arguments=arguments,
            check=check,
            location=None,
            pipe_input=pipe_input,
        )

    def is_available(self):
        # try:
        #     result = self.execute(
        #         command=EE.OMEGA, arguments=[EE.OMEGA_HELP], check=True
        #     )

        #     if EE.OMEGA_HELP_IDENTIFICATION_STRING in result.stderr:
        #         return True
        #     return False
        # except Exception as e:
        #     return False
        pass
