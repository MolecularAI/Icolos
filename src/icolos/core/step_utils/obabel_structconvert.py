from icolos.utils.enums.program_parameters import OpenBabelEnum
from icolos.utils.execute_external.openbabel import OpenBabelExecutor
from icolos.loggers.steplogger import StepLogger
from icolos.utils.enums.program_parameters import SchrodingerExecutablesEnum
from icolos.utils.execute_external.schrodinger import SchrodingerExecutor
from icolos.utils.general.icolos_exceptions import StepFailed
from icolos.utils.enums.logging_enums import LoggingConfigEnum

_LE = LoggingConfigEnum()
_OE = OpenBabelEnum()


class OBabelStructConvert:
    """
    Utility for converting structure files with Schrodinger's StructConvert
    """

    def __init__(
        self,
    ) -> None:
        self._logger = StepLogger()

        self.executor = OpenBabelExecutor()
        if not self.executor.is_available():
            # this shouldn't fail, obabel is in the env
            raise StepFailed("Cannot initialize sdconvert backend - abort.")
        self._logger.log(f"Checked obabel availability - valid.", _LE.DEBUG)

    def execute(self, arguments: list):
        execution_result = self.executor.execute(
            command=_OE.OBABEL, arguments=arguments, check=True
        )
        if execution_result.returncode != 0:
            self._logger.log(
                f"Could not execute sdconvert (returncode != 0) with error: {execution_result.stderr}.",
                _LE.ERROR,
            )

    def convert(self, input_file: str, output_file: str):
        arguments = [
            input_file,
            output_file,
        ]
        self.execute(arguments=arguments)

    def sdf2pdb(self, sdf_file: str, pdb_file: str):
        assert sdf_file.endswith(".sdf")
        assert pdb_file.endswith(".pdb")

        arguments = ["-isdf", sdf_file, "-O", pdb_file]
        self.execute(arguments=arguments)
