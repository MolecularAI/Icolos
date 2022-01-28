from icolos.loggers.steplogger import StepLogger
from icolos.utils.execute_external.sdconvert import SDConvertExecutor

from icolos.utils.enums.program_parameters import SchrodingerExecutablesEnum
from icolos.utils.enums.logging_enums import LoggingConfigEnum
from icolos.utils.general.icolos_exceptions import StepFailed

_SEE = SchrodingerExecutablesEnum()
_LE = LoggingConfigEnum()


class SDConvertUtil:
    def __init__(self, prefix_execution: str = None, binary_location: str = None):
        self._logger = StepLogger()

        # initialize and check executor
        self.executor = SDConvertExecutor(
            prefix_execution=prefix_execution, binary_location=binary_location
        )
        if not self.executor.is_available():
            raise StepFailed("Cannot initialize sdconvert backend - abort.")
        self._logger.log(f"Checked sdconvert availability - valid.", _LE.DEBUG)

    def execute(self, arguments: list):
        execution_result = self.executor.execute(
            command=_SEE.SDCONVERT, arguments=arguments, check=True
        )
        if execution_result.returncode != 0:
            self._logger.log(
                f"Could not execute sdconvert (returncode != 0) with error: {execution_result.stderr}.",
                _LE.ERROR,
            )

    def mae2sdf(self, mae_file: str, sdf_file: str):
        arguments = [
            "".join([_SEE.SDCONVERT_I, _SEE.SDCONVERT_FORMAT_MAE]),
            mae_file,
            "".join([_SEE.SDCONVERT_O, _SEE.SDCONVERT_FORMAT_SD]),
            sdf_file,
        ]
        self.execute(arguments=arguments)

    def sdf2mae(self, sdf_file: str, mae_file: str):
        arguments = [
            "".join([_SEE.SDCONVERT_I, _SEE.SDCONVERT_FORMAT_SD]),
            sdf_file,
            "".join([_SEE.SDCONVERT_O, _SEE.SDCONVERT_FORMAT_MAE]),
            mae_file,
        ]
        self.execute(arguments=arguments)

    def pdb2mae(self, pdb_file: str, mae_file: str):
        arguments = [
            "".join([_SEE.SDCONVERT_I, _SEE.SDCONVERT_FORMAT_SD]),
            pdb_file,
            "".join([_SEE.SDCONVERT_O, _SEE.SDCONVERT_FORMAT_MAE]),
            mae_file,
        ]
        self.execute(arguments=arguments)

    def sdf2pdb(self, sdf_file: str, pdb_file: str):
        arguments = [
            "".join([_SEE.SDCONVERT_I, _SEE.SDCONVERT_FORMAT_SD]),
            pdb_file,
            "".join([_SEE.SDCONVERT_O, _SEE.SDCONVERT_FORMAT_PDB]),
            sdf_file,
        ]
        self.execute(arguments=arguments)
