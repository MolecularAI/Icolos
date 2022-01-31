from icolos.loggers.steplogger import StepLogger
from icolos.utils.enums.program_parameters import SchrodingerExecutablesEnum
from icolos.utils.execute_external.schrodinger import SchrodingerExecutor
from icolos.utils.general.icolos_exceptions import StepFailed
from icolos.utils.enums.logging_enums import LoggingConfigEnum

_LE = LoggingConfigEnum()
_SEE = SchrodingerExecutablesEnum()


class StructConvert:
    """
    Utility for converting structure files with Schrodinger's StructConvert
    """

    def __init__(self, prefix_execution: str, binary_location: str = None) -> None:
        self._logger = StepLogger()

        self.executor = SchrodingerExecutor(
            binary_location=binary_location, prefix_execution=prefix_execution
        )
        if not self.executor.is_available():
            raise StepFailed("Cannot initialize sdconvert backend - abort.")
        self._logger.log(f"Checked sdconvert availability - valid.", _LE.DEBUG)

    def execute(self, arguments: list):
        execution_result = self.executor.execute(
            command=_SEE.STRUCTCONVERT, arguments=arguments, check=True
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

    def pdb2mae(self, pdb_file: str, mae_file: str):
        # new schrodinger does not check this, needs to be done here
        assert pdb_file.endswith(".pdb")
        assert mae_file.endswith(".mae")
        arguments = [
            pdb_file,
            mae_file,
        ]
        self.execute(arguments=arguments)

    def sdf2pdb(self, sdf_file: str, pdb_file: str):
        assert sdf_file.endswith(".sdf")
        assert pdb_file.endswith(".pdb")
        arguments = [
            sdf_file,
            pdb_file,
        ]
        self.execute(arguments=arguments)

    def mae2pdb(self, mae_file: str, pdb_file: str):
        assert mae_file.endswith(".mae")
        assert pdb_file.endswith(".pdb")
        arguments = [
            mae_file,
            pdb_file,
        ]
        self.execute(arguments=arguments)
