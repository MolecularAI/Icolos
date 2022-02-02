from icolos.core.workflow_steps.schrodinger.fep_base import StepFEPBase
from icolos.utils.enums.step_enums import StepBaseEnum, StepFepPlusEnum
from icolos.utils.enums.program_parameters import FepPlusEnum

from pydantic import BaseModel

_FE = FepPlusEnum()
_SFE = StepFepPlusEnum()
_SBE = StepBaseEnum


class StepFepPlusAnalysis(StepFEPBase, BaseModel):
    """
    Standalone class to analyse data from a previous fep job
    """

    def __init__(self, **data):
        super().__init__(**data)

    def execute(self):
        """
        Analyses the map produced from an FEP run
        """
        tmp_dir = self._make_tmpdir()
        self.data.generic.write_out_all_files(tmp_dir)
        self._extract_log_file_data(tmp_dir)
        self._remove_temporary(tmp_dir)
