from typing import Dict
from pydantic import BaseModel
import os
import sys

from icolos.core.composite_agents.workflow import WorkFlow
from icolos.utils.enums.logging_enums import LoggingConfigEnum

from icolos.utils.entry_point_functions.logging_helper_functions import (
    initialize_logging,
)
from icolos.utils.entry_point_functions.parsing_functions import (
    get_runtime_global_variables,
    add_global,
)
from icolos.utils.general.files_paths import attach_root_path
from icolos.utils.enums.composite_agents_enums import WorkflowEnum

_WE = WorkflowEnum()
_LE = LoggingConfigEnum()


class IcolosWorkflow(BaseModel):
    """
    Alternative programmatic entrypoint for the Icolos workflow
    """

    config: Dict = None
    workflow: WorkFlow = None
    logging: str = None
    global_vars: Dict = None

    def __init__(self, config, global_vars: Dict = None) -> None:
        super().__init__(**config)

        self.config = self._parse_global_vars(config, global_vars)
        self.workflow = WorkFlow(**config[_WE.WORKFLOW])
        # tutorial settings logs everything to stdout as well as the file
        self.logging = "tutorial"

    def _initialize_logging(self):
        log_conf = attach_root_path(_LE.PATH_CONFIG_TUTORIAL)
        logger = initialize_logging(log_conf_path=log_conf, workflow_conf=self.config)
        return logger

    def _parse_global_vars(self, config, global_vars):
        # substitute global vars throughout the config file, return modified config

        if global_vars is not None:
            config = add_global(config, global_vars, _WE.GLOBAL_VARIABLES)
        config = add_global(
            config,
            get_runtime_global_variables(
                os.path.join(os.getcwd(), "config.json"), os.path.realpath(__file__)
            ),
            _WE.GLOBAL_VARIABLES,
        )
        return config

    def execute(self):
        self._initialize_logging()
        self.workflow.initialize()
        self.workflow.execute()

        sys.exit(0)
