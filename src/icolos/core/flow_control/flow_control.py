from typing import List
from pydantic import BaseModel, PrivateAttr
from icolos.core.workflow_steps.step import StepSettingsParameters
from icolos.core.workflow_steps.step import StepBase
from icolos.loggers.steplogger import StepLogger
from icolos.core.workflow_steps.step import (
    StepData,
    StepInputParameters,
    StepWriteoutParameters,
    StepExecutionParameters,
)
from icolos.utils.enums.step_enums import StepBaseEnum
from icolos.utils.enums.step_initialization_enum import StepInitializationEnum
from icolos.utils.general.convenience_functions import nested_get

_SIE = StepInitializationEnum()


class BaseStepConfig(BaseModel):
    """
    Minimal template class for the base config, without the unnecessary stuff that StepBase requires
    """

    step_id: str = None
    work_dir: str = None
    type: str = None
    data: StepData = StepData()
    input: StepInputParameters = StepInputParameters()
    writeout: List[StepWriteoutParameters] = []
    execution: StepExecutionParameters = StepExecutionParameters()
    settings: StepSettingsParameters = StepSettingsParameters()

    def _as_dict(self):
        return {
            "step_id": self.step_id,
            "type": self.type,
            "execution": self.execution,
            "settings": self.settings,
            "work_dir": self.work_dir,
            "data": self.data,
            "input": self.input,
            "writeout": self.writeout,
        }


class FlowControlBase(BaseModel):
    # List of steps to be iterated over, each set needs their inputs chained together
    base_config: List[BaseStepConfig] = None
    initialized_steps: List[StepBase] = None
    _logger = PrivateAttr()

    def __init__(self, **data) -> None:
        super().__init__(**data)
        self._logger = StepLogger()

    def _initialize_step_from_dict(self, step_conf: dict):
        # Require a separate initialisation method to avoid circular import
        _STE = StepBaseEnum

        step_type = nested_get(step_conf, _STE.STEP_TYPE, default=None)
        step_type = None if step_type is None else step_type.upper()
        if step_type in _SIE.STEP_INIT_DICT.keys():
            return _SIE.STEP_INIT_DICT[step_type](**step_conf)
        else:
            raise ValueError(
                f"Backend for step {nested_get(step_conf, _STE.STEPID, '')} unknown."
            )
