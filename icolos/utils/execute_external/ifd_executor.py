from icolos.utils.execute_external.execute import ExecutorBase
from icolos.utils.enums.program_parameters import InducedFitEnum

_IFE = InducedFitEnum()


class IFDExecutor(ExecutorBase):
    def __init__(self, prefix_execution=None, binary_location=None) -> None:
        super().__init__(
            prefix_execution=prefix_execution, binary_location=binary_location
        )

    def execute(
        self, command: str, arguments: list, check=True, location=None, pipe_input=None
    ):

        if command not in [_IFE.IFD_EXEC]:
            raise AssertionError(
                "Commmand must be recognised in the internal dictionary"
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
                command=_IFE.IFD_EXEC, arguments=[_IFE.IFD_HELP], check=True
            )

            if _IFE.IFD_HELP_ID in result.stdout:
                return True
            return False
        except Exception as e:
            return False
