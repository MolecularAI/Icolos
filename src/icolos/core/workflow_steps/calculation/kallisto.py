from typing import List
from pydantic import BaseModel

from icolos.core.workflow_steps.calculation.base import StepCalculationBase
from icolos.utils.enums.program_parameters import KallistoEnum
from icolos.utils.enums.step_enums import StepKallistoEnum
from icolos.core.workflow_steps.step import _LE
from icolos.utils.execute_external.kallisto import KallistoExecutor

_SKE = StepKallistoEnum()
_KE = KallistoEnum()


class KallistoAdditional(BaseModel):
    features: List[str] = [_KE.ALP, _KE.BONDS]        # list of features to be obtained


class StepKallisto(StepCalculationBase, BaseModel):

    kallisto_additional: KallistoAdditional = None

    def __init__(self, **data):
        super().__init__(**data)

        self._initialize_backend(executor=KallistoExecutor)
        self._check_backend_availability()

        # initialize the additional settings
        self.kallisto_additional = KallistoAdditional(**self.settings.additional)

    def execute(self):
        """tmp_dir = self._make_tmpdir()
        self._write_panther_config_file(tmp_dir)
        self._execute_backend(tmp_dir)
        self._logger.log("Executed PANTHER and obtained negative image.", _LE.INFO)
        self._logger.log(
            f"Calculated negative image for configuration file in {tmp_dir}.", _LE.DEBUG
        )
        self._parse_panther_output(tmp_dir)
        self._remove_temporary(tmp_dir)"""
        print(self.kallisto_additional)
