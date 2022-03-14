from typing import List

from icolos.loggers.steplogger import StepLogger
from icolos.utils.execute_external.openbabel import OpenBabelExecutor
from icolos.utils.execute_external.structcat import StructcatExecutor

from icolos.utils.enums.program_parameters import (
    OpenBabelEnum,
    SchrodingerExecutablesEnum,
)
from icolos.utils.enums.logging_enums import LoggingConfigEnum
from icolos.utils.general.icolos_exceptions import StepFailed

_SEE = SchrodingerExecutablesEnum()
_LE = LoggingConfigEnum()

_OE = OpenBabelEnum()


class StructcatUtil:
    def __init__(
        self,
        prefix_execution: str = None,
        binary_location: str = None,
        backend: str = _OE.OBABEL,
    ):
        self._logger = StepLogger()
        self._backend = backend
        # initialize and check executor
        if self._backend == _SEE.STRUCTCAT:
            self.executor = StructcatExecutor(
                prefix_execution=prefix_execution, binary_location=binary_location
            )
        elif self._backend == _OE.OBABEL:
            self.executor = OpenBabelExecutor()

        if not self.executor.is_available():
            raise StepFailed("Cannot initialize structcat backend - abort.")

    def concatenate(
        self,
        input_files: List[str],
        output_file: str,
        location: str = None,
    ):
        if self._backend == _SEE.STRUCTCAT:
            arguments = []
            for input_file in input_files:
                arguments = arguments + [
                    _SEE.STRUCTCAT_I,
                    input_file,
                ]
            arguments = arguments + [
                _SEE.STRUCTCAT_O,
                output_file,
            ]
            self.executor.execute(
                command=_SEE.STRUCTCAT, arguments=arguments, check=True
            )

        elif self._backend == _OE.OBABEL:
            arguments = input_files
            arguments.append("-O")
            arguments.append(output_file)
            self.executor.execute(
                command=_OE.OBABEL, arguments=arguments, check=True, location=location
            )
