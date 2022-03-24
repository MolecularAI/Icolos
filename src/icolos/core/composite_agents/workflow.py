from typing import Dict, List

from pydantic import BaseModel, PrivateAttr
from icolos.core.containers.gmx_state import GromacsState
from icolos.core.containers.perturbation_map import PerturbationMap
from icolos.core.flow_control.flow_control import FlowControlBase
from icolos.core.workflow_steps.step import StepBase
from icolos.core.composite_agents.base_agent import BaseAgent, AgentHeaderParameters
from icolos.utils.enums.step_enums import StepBaseEnum

from icolos.utils.general.icolos_exceptions import get_exception_message

from icolos.utils.enums.logging_enums import LoggingConfigEnum

_LE = LoggingConfigEnum()
_SBE = StepBaseEnum


class WorkflowHeaderParameters(AgentHeaderParameters, BaseModel):
    pass


class WorkflowData(BaseModel):
    work_dir: str = None
    perturbation_map: PerturbationMap = None


class WorkFlow(BaseAgent, BaseModel):
    """Class to hold the whole logic for a workflow."""

    steps: List[Dict] = []
    header: WorkflowHeaderParameters = WorkflowHeaderParameters()
    workflow_data: WorkflowData = WorkflowData()

    class Config:
        underscore_attrs_are_private = True

    _logger = PrivateAttr()
    _initialized_steps = PrivateAttr()

    def __init__(self, **data):
        super().__init__(**data)
        self._initialized_steps = []

    def initialize(self):
        # this feels hacky but it fixes a circular import and allows a step in the workflow
        # to create other workflows
        from icolos.core.steps_utils import initialize_step_from_dict

        super().initialize()
        self._initialized_steps = []
        for step_conf in self.steps:
            step_conf = self._update_global_variables(conf=step_conf)
            step = initialize_step_from_dict(step_conf=step_conf)
            if isinstance(step, StepBase):
                # we have a normal step, no flow control wrapping
                step.set_workflow_object(self)
                self._initialized_steps.append(step)
            elif isinstance(step, FlowControlBase):
                self._initialized_steps.append(step.dispatcher)
                step.dispatcher.set_workflow_object(self)
        self._logger.log(
            f"Initialized {len(self._initialized_steps)} steps in workflow {self.header.id}.",
            _LE.DEBUG,
        )

    def execute(self):
        for step in self._initialized_steps:
            step.generate_input()
            self._logger.log(f"Starting execution of step: {step.step_id}", _LE.INFO)
            step.execute()
            self._logger.log(
                f"Processing write-out blocks for {step.step_id}.", _LE.DEBUG
            )
            step.process_write_out()
        self._logger.log(
            f"Execution of {len(self._initialized_steps)} steps completed.", _LE.INFO
        )

    def is_valid(self) -> bool:
        try:
            for step in self._initialized_steps:
                step.validate()
        except Exception as e:
            self._logger.log(
                f'During step validation, "WorkFlow" returned the following exception: {get_exception_message(e)}.',
                _LE.WARNING,
            )
            return False
        return True

    def add_step(self, step: StepBase):
        self._initialized_steps.append(step)

    def get_steps(self) -> list:
        return self._initialized_steps

    def find_step_by_step_id(self, step_id: str):
        for step in self._initialized_steps:
            if step.step_id == step_id:
                return step
            elif step.step_id == _SBE.STEP_DISPATCHER:
                pass

        raise IndexError(f"Could not find step with step_id {step_id} in workflow.")

    def __iter__(self):
        return iter(self.steps)

    def __repr__(self):
        return "<Icolos workflow: id=%s, description=%s, number steps: %s>" % (
            self.get_id(),
            self.get_description(),
            len(self),
        )

    def set_perturbation_map(self, p_map: PerturbationMap) -> None:
        self.workflow_data.perturbation_map = p_map

    def get_perturbation_map(self) -> PerturbationMap:
        return self.workflow_data.perturbation_map

    def __str__(self):
        return self.__repr__()

    def __getitem__(self, key: int):
        return self._initialized_steps[key]

    def __len__(self) -> int:
        return len(self._initialized_steps)
