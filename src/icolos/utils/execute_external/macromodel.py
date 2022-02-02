from icolos.utils.enums.program_parameters import (
    MacromodelEnum,
    SchrodingerExecutablesEnum,
)
from icolos.utils.execute_external.execute import ExecutorBase

SEE = SchrodingerExecutablesEnum()
EE = MacromodelEnum()


class MacromodelExecutor(ExecutorBase):
    """For the execution of the "macromodel" binary."""

    def __init__(self, prefix_execution=None, binary_location=None):
        super().__init__(
            prefix_execution=prefix_execution, binary_location=binary_location
        )

    def execute(
        self, command: str, arguments: list, check=True, location=None, pipe_input=None
    ):
        # check, whether a proper executable is provided
        if command not in [EE.MACROMODEL]:
            raise ValueError(
                "Parameter command must be an dictionary of the internal Macromodel executable list."
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
                command=EE.MACROMODEL, arguments=[EE.MACROMODEL_HELP], check=True
            )

            if EE.MACROMODEL_HELP_IDENTIFICATION_STRING in result.stdout:
                return True
            return False
        except Exception as e:
            return False
