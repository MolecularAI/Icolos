from pydantic import BaseModel
from icolos.core.workflow_steps.step import StepBase
from icolos.utils.enums.step_enums import StepBaseEnum

_SBE = StepBaseEnum


class StepIOBase(StepBase, BaseModel):
    def __init__(self, **data):
        super().__init__(**data)
