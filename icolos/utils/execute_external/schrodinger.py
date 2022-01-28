from icolos.utils.enums.program_parameters import SchrodingerExecutablesEnum
from icolos.utils.execute_external.execute import ExecutorBase

_SEE = SchrodingerExecutablesEnum()


class SchrodingerExecutor(ExecutorBase):
    """For the execution of Schrodinger's support entry points"""

    def __init__(self, prefix_execution=None, binary_location=None):
        super().__init__(
            prefix_execution=prefix_execution, binary_location=binary_location
        )

    def execute(
        self, command: str, arguments: list, check=True, location=None, pipe_input=None
    ):
        # check, whether a proper executable is provided; update the calls to "$SCHRODINGER/XYZ"
        if command == _SEE.STRUCTCAT:
            command = _SEE.STRUCTCAT_CALL
        elif command == _SEE.SDCONVERT:
            command = _SEE.SDCONVERT_CALL
        elif command == _SEE.STRUCT_SPLIT:
            command = _SEE.STRUCT_SPLIT_CALL
        elif command == _SEE.STRUCTCONVERT:
            command = _SEE.STRUCTCONVERT_CALL
        elif command == _SEE.FMP_STATS:
            command = _SEE.FMP_STATS_CALL
        elif command == _SEE.PREPWIZARD:
            command = _SEE.PREPWIZARD_CALL
        elif command == _SEE.MULTISIM_EXEC:
            command = _SEE.MULTISIM_EXEC
        else:
            raise ValueError(
                "Parameter command must be an dictionary of the internal Schrodinger entry point list."
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
                command=_SEE.STRUCTCAT, arguments=[_SEE.STRUCTCAT_HELP], check=True
            )

            if _SEE.STRUCTCAT_HELP_IDENTIFICATION_STRING in result.stdout:
                return True
            return False
        except Exception as e:

            print(e)
            return False
