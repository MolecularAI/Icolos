from icolos.core.flow_control.iterator import StepIterator
from icolos.core.workflow_steps.active_learning.active_learning_reinvent import (
    StepActiveLearning,
)
from icolos.utils.enums.step_enums import StepBaseEnum

_SBE = StepBaseEnum


class FlowControlInitializationEnum:
    # These steps are responsible for initializing other steps as part of their execution
    # Keep these separate to the main pool of steps to avoid circular imports

    FLOW_CONTROL_INIT_DICT = {
        _SBE.STEP_ACTIVE_LEARNING: StepActiveLearning,
        _SBE.STEP_ITERATOR: StepIterator,
    }
