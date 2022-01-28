from icolos.utils.execute_external.execute import ExecutorBase
from icolos.utils.enums.program_parameters import (
    FepPlusEnum,
    SchrodingerExecutablesEnum,
)

FE = FepPlusEnum()
SEE = SchrodingerExecutablesEnum()


class FepPlusExecutor(ExecutorBase):
    def __init__(self, prefix_execution=None, binary_location=None):
        super().__init__(
            prefix_execution=prefix_execution, binary_location=binary_location
        )

    def execute(
        self, command: str, arguments: list, check=True, location=None, pipe_input=None
    ):
        if command not in [
            FE.FEP_MAPPER,
            FE.FEP_EXECUTOR,
            FE.JSC_LIST,
            FE.JSC_TAIL_FILE,
        ]:
            raise ValueError(
                "Execution command must be recognised by the executable's enum"
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
                command=FE.FEP_MAPPER, arguments=[FE.FEP_HELP], check=True
            )
            if FE.FEP_MAPPER_HELP_SUCCESS_STRING in result.stdout:
                return True
            return False
        except Exception as e:
            print(str(e))
            return False
