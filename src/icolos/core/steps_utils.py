from icolos.core.workflow_steps.step import StepBase
from icolos.utils.enums.step_enums import StepBaseEnum
from icolos.utils.general.convenience_functions import nested_get
from icolos.utils.enums.step_initialization_enum import StepInitializationEnum
from icolos.utils.enums.flow_control_enums import FlowControlInitializationEnum

_IE = StepInitializationEnum()
_FCE = FlowControlInitializationEnum()


def initialize_step_from_dict(step_conf: dict) -> StepBase:
    _STE = StepBaseEnum
    step_type = nested_get(step_conf, _STE.STEP_TYPE, default=None)
    step_type = None if step_type is None else step_type.upper()
    if step_type in _IE.STEP_INIT_DICT.keys():
        return _IE.STEP_INIT_DICT[step_type](**step_conf)
    elif step_type in _FCE.FLOW_CONTROL_INIT_DICT.keys():
        return _FCE.FLOW_CONTROL_INIT_DICT[step_type](**step_conf)
    else:
        raise ValueError(
            f"Backend for step {nested_get(step_conf, _STE.STEPID, '')} unknown."
        )
