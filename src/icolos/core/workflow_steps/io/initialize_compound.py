from pydantic import BaseModel

from icolos.utils.general.icolos_exceptions import StepFailed
from icolos.core.workflow_steps.io.base import StepIOBase
from icolos.core.workflow_steps.step import _LE


class StepInitializeCompound(StepIOBase, BaseModel):
    def __init__(self, **data):
        super().__init__(**data)

    def execute(self):
        if len(self.data.compounds) == 0:
            raise StepFailed(
                "Compound initialization step failed - no Compound objects generated."
            )
        self._logger.log(
            f"Step {self.get_step_id()} initialized {len(self.get_compounds())} compounds.",
            _LE.INFO,
        )
