from icolos.utils.enums.program_parameters import AutoDockVinaEnum
from icolos.utils.execute_external.execute import ExecutorBase

_EE = AutoDockVinaEnum()


class AutoDockVinaExecutor(ExecutorBase):
    """For the execution of AutoDock Vina 1.2.0."""

    def __init__(self, prefix_execution=None, binary_location=None):
        super().__init__(
            prefix_execution=prefix_execution, binary_location=binary_location
        )

    def execute(
        self, command: str, arguments: list, check=True, location=None, pipe_input=None
    ):
        # check, whether a proper executable is provided
        if command not in [_EE.VINA_CALL]:
            raise ValueError(
                "Parameter command must be in dictionary of the internal AutoDock Vina executable list."
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
            result = self.execute(
                command=_EE.VINA_CALL, arguments=[_EE.VINA_HELP], check=True
            )
            if result.returncode == 0:
                return True
            return False
        except Exception as e:
            return False
