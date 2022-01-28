from pydantic import BaseModel
from icolos.core.workflow_steps.step import StepBase


class StepConfgenBase(StepBase, BaseModel):
    def __init__(self, **data):
        super().__init__(**data)
