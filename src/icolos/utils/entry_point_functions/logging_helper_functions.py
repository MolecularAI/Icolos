import json
import logging.config as logging_config
from icolos.loggers.entrypoint_logger import EntryPointLogger
from icolos.utils.enums.composite_agents_enums import WorkflowEnum
from icolos.utils.enums.logging_enums import LoggingConfigEnum

from icolos.utils.general.convenience_functions import *

_WE = WorkflowEnum()
_LE = LoggingConfigEnum()


def initialize_logging(log_conf_path: str, workflow_conf: dict) -> EntryPointLogger:
    with open(log_conf_path, "r") as f:
        log_conf_dict = json.load(f)
    header = nested_get(workflow_conf, [_WE.WORKFLOW, _WE.HEADER], default={})
    if in_keys(header, [_WE.LOGGING, _WE.LOGGING_LOGFILE]):
        try:
            log_conf_dict["handlers"]["file_handler"]["filename"] = nested_get(
                header, [_WE.LOGGING, _WE.LOGGING_LOGFILE], None
            )
            log_conf_dict["handlers"]["file_handler_blank"]["filename"] = nested_get(
                header, [_WE.LOGGING, _WE.LOGGING_LOGFILE], None
            )
        except KeyError:
            pass
    logging_config.dictConfig(log_conf_dict)
    logger = EntryPointLogger()
    return logger
