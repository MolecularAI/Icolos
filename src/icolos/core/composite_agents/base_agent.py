import os
from abc import abstractmethod
from copy import deepcopy
from typing import Dict, List

from pydantic import BaseModel, PrivateAttr

from icolos.loggers.agentlogger import AgentLogger

from icolos.utils.enums.logging_enums import LoggingConfigEnum
from icolos.utils.enums.composite_agents_enums import WorkflowEnum

_WE = WorkflowEnum()
_LE = LoggingConfigEnum()


class AgentEnvironmentParameters(BaseModel):
    class WorkflowExportParameters(BaseModel):
        key: str
        value: str

    export: List[WorkflowExportParameters] = []


class AgentHeaderParametersSettings(BaseModel):
    remove_temporary_files: bool = True
    single_directory: bool = False


class AgentHeaderParameters(BaseModel):
    class AgentLoggingParameters(BaseModel):
        logfile: str = None

    id: str = None
    description: str = None
    logging: AgentLoggingParameters = AgentLoggingParameters()
    environment: AgentEnvironmentParameters = None
    global_variables: Dict = None
    global_settings: AgentHeaderParametersSettings = AgentHeaderParametersSettings()


class BaseAgent(BaseModel):

    # should also work without parsing the base specification here, but then IDEs will not pick up stuff below
    header: AgentHeaderParameters = AgentHeaderParameters()

    class Config:
        underscore_attrs_are_private = True

    _logger = PrivateAttr()

    def __init__(self, **data):
        super().__init__(**data)
        self._logger = AgentLogger()

    def _export_env_variables(self):
        for var in self.header.environment.export:
            key = str(var.key)
            value = os.path.expandvars(str(var.value))
            os.environ[key] = value
            self._logger.log(f"Exported variable {key} with value {value}.", _LE.DEBUG)

    def initialize(self):
        self._export_env_variables()

    def _nested_update(self, inp, pattern: str, replacement: str):
        if isinstance(inp, dict):
            items = inp.items()
        elif isinstance(inp, (list, tuple)):
            items = enumerate(inp)
        elif isinstance(inp, str):
            return inp.replace(pattern, replacement)
        else:
            return inp

        for key, value in items:
            inp[key] = self._nested_update(value, pattern, replacement)
        return inp

    def _update_global_variables(self, conf: dict) -> dict:
        conf = deepcopy(conf)
        if self.header.global_variables is not None:
            for key, value in self.header.global_variables.items():
                pattern = "{" + key + "}"
                self._nested_update(inp=conf, pattern=pattern, replacement=value)
                self._logger.log(
                    f"Updated global variable {key} with value {value}.", _LE.DEBUG
                )
        return conf

    @abstractmethod
    def execute(self):
        raise NotImplementedError

    def is_valid(self) -> bool:
        raise NotImplementedError

    def set_id(self, id: str):
        self.header.id = id

    def get_id(self) -> str:
        return self.header.id

    def set_description(self, description: str):
        self.header.description = description

    def get_description(self) -> str:
        return self.header.description
