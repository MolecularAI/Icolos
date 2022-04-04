import os
from typing import Tuple

from icolos.loggers.base_logger import BaseLogger

from icolos.utils.enums.composite_agents_enums import WorkflowEnum
from icolos.utils.enums.entry_points import ExecutorEnum
from icolos.utils.enums.logging_enums import LoggingConfigEnum
from icolos.utils.general.files_paths import move_up_directory

_WE = WorkflowEnum()
_LE = LoggingConfigEnum()
_EE = ExecutorEnum()


def parse_global(g_input, logger: BaseLogger) -> dict:
    if g_input is not None:
        if not isinstance(g_input, list):
            g_input = [g_input]
        g_vars = {}
        for new_var in g_input:
            parts = new_var.split(":")
            if len(parts) != 2:
                logger.log(
                    f"Ignoring global input {new_var} set by command-line, as they must have one key and one value, separated by ':'.",
                    _LE.WARNING,
                )
                continue
            g_vars[parts[0]] = parts[1]
            logger.log(
                f'Parsed global input "{parts[0]}" (value: "{parts[1]}").', _LE.DEBUG
            )
        return g_vars
    else:
        return {}


def add_global(configuration: dict, g_vars: dict, field: str) -> dict:
    """This function adds (and overwrites) values for global settings and variables. Parameter "field" selects,
    which key is to be used in the header region."""
    header = configuration[_WE.WORKFLOW][_WE.HEADER]
    if field not in header.keys():
        header[field] = {}
    for key, value in g_vars.items():
        header[field][key] = value
    return configuration


def get_runtime_global_variables(args_conf: str, entry_point_path: str) -> dict:
    return {  # current workdir
        _EE.RUNTIME_GLOBAL_VARIABLE_WORKDIR: os.getcwd(),
        # directory where the entry point lies
        _EE.RUNTIME_GLOBAL_VARIABLE_ENTRYPOINTDIR: os.path.dirname(entry_point_path),
        # directory where the JSON lies
        _EE.RUNTIME_GLOBAL_VARIABLE_CONFIGDIR: os.path.dirname(
            os.path.abspath(args_conf)
        ),
        # top level of Icolos icolos package
        _EE.RUNTIME_GLOBAL_VARIABLE_PACKAGEDIR: move_up_directory(
            os.path.dirname(entry_point_path), 3
        ),
    }


def parse_header(conf: dict, args, entry_point_path: str, logger: BaseLogger) -> dict:
    # parse global variables from command-line
    global_vars_CLI = parse_global(g_input=args.global_variables, logger=logger)
    conf = add_global(conf, global_vars_CLI, _WE.GLOBAL_VARIABLES)

    # add run-specified global variables (the current directory, the JSONs directory, ...)
    conf = add_global(
        conf,
        get_runtime_global_variables(args.conf, entry_point_path),
        _WE.GLOBAL_VARIABLES,
    )

    # update global settings; if they are not supported, pydantic will complain later on
    # TODO: at the moment this implementation ignores stuff that is not understood (e.g. when a typo occurs); this should fail
    global_settings_CLI = parse_global(g_input=args.global_settings, logger=logger)
    conf = add_global(conf, global_settings_CLI, _WE.GLOBAL_SETTINGS)
    return conf


def get_version_number() -> str:
    try:
        # this requires python >= 3.8
        from importlib import metadata

        return metadata.version("icolos")
    except:
        return None


def get_config_version_number(conf: dict) -> str:
    try:
        return str(conf[_WE.WORKFLOW][_WE.HEADER][_WE.VERSION])
    except:
        return None


def log_version_number(logger: BaseLogger):
    version = get_version_number()
    if version is None:
        logger.log(f"Could not obtain Icolos version.", _LE.WARNING)
        version = "unknown"
    logger.log(f"Icolos version {version} initialized.", _LE.INFO)


def version_match(conf: dict) -> bool:
    version_config = get_config_version_number(conf)
    version_installation = get_version_number()
    if (
        version_config is None
        or version_installation is None
        or version_config != version_installation
    ):
        return False
    return True


def get_versions_as_strings(conf: dict) -> Tuple[str, str]:
    version_config = get_config_version_number(conf)
    version_installation = get_version_number()
    version_config = "unknown" if version_config is None else version_config
    version_installation = (
        "unknown" if version_installation is None else version_installation
    )
    return version_config, version_installation
