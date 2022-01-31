from pydantic import BaseModel, PrivateAttr

from icolos.core.composite_agents.base_agent import BaseAgent, AgentHeaderParameters

from icolos.utils.enums.logging_enums import LoggingConfigEnum
from icolos.utils.enums.composite_agents_enums import SchedulerEnum

_SE = SchedulerEnum()
_LE = LoggingConfigEnum()


class SchedulerHeaderParameters(AgentHeaderParameters, BaseModel):
    pass


class Scheduler(BaseAgent, BaseModel):
    """Class to hold the whole logic for scheduling sub-jobs."""

    header: SchedulerHeaderParameters = SchedulerHeaderParameters()

    class Config:
        underscore_attrs_are_private = True

    _logger = PrivateAttr()
    _initialized_steps = PrivateAttr()

    def __init__(self, **data):
        super().__init__(**data)

    def initialize(self):
        super().initialize()

    def execute(self):
        # TODO: implement
        pass

    def _action_prepare(self):
        pass

    def _action_run(self):
        pass

    def is_valid(self) -> bool:
        # TODO: implement
        pass

    def __repr__(self):
        return "<Icolos scheduler: id=%s, description=%s>" % (
            self.get_id(),
            self.get_description(),
        )

    def __str__(self):
        return self.__repr__()
