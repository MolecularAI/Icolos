from icolos.core.flow_control.iterator import StepIterator
from icolos.utils.enums.step_enums import StepBaseEnum

_SBE = StepBaseEnum


class FlowControlInitializationEnum:
    # These steps are responsible for initializing other steps as part of their execution
    # Keep these separate to the main pool of steps to avoid circular imports

    FLOW_CONTROL_INIT_DICT = {
        _SBE.STEP_ITERATOR: StepIterator,
    }
