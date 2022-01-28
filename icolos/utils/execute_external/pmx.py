from icolos.utils.enums.program_parameters import PMXEnum
from icolos.utils.execute_external.execute import ExecutorBase

_PE = PMXEnum()


class PMXExecutor(ExecutorBase):
    def __init__(self, prefix_execution=None, binary_location=None):
        super().__init__(
            prefix_execution=prefix_execution, binary_location=binary_location
        )

    def execute(
        self, command: str, arguments: list, check=True, location=None, pipe_input=None
    ):
        if command not in [
            _PE.ABFE,
            _PE.ANALYSE,
            _PE.ATOMMAPPING,
            _PE.DOUBLEBOX,
            _PE.GENLIB,
            _PE.GENTOP,
            _PE.LIGANDHYBRID,
            _PE.MUTATE,
            _PE.BOX_WATER_IONS,
            _PE.PREPARE_SIMULATIONS,
            _PE.PREPARE_TRANSITIONS,
            _PE.RUN_ANALYSIS,
            _PE.RUN_SIMULATIONS,
            _PE.ASSEMBLE_SYSTEMS,
        ]:
            raise ValueError(
                "Command must be present in internal list of PMX executables."
            )

        # handle for dealing with programs that want interactive input
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
                command=_PE.ANALYSE, arguments=[_PE.ANALYSE_HELP], check=False
            )
            if _PE.ANALYSE_HELP_SUCCESS_STRING in result.stdout:
                return True
            return False
        except Exception as e:
            return False
