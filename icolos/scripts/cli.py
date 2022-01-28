#!/usr/bin/env python

import os
import sys
import json
import argparse

from icolos.core.composite_agents.workflow import WorkFlow

from icolos.loggers.entrypoint_logger import EntryPointLogger

from icolos.utils.enums.composite_agents_enums import WorkflowEnum
from icolos.utils.enums.logging_enums import LoggingConfigEnum
from icolos.utils.enums.entry_points import ExecutorEnum

from icolos.utils.entry_point_functions.logging_helper_functions import (
    initialize_logging,
)
from icolos.utils.entry_point_functions.parsing_functions import parse_header
from icolos.utils.general.files_paths import attach_root_path


class IcolosCLI:
    def __init__(self) -> None:
        # enums
        _LE = LoggingConfigEnum()
        _EE = ExecutorEnum()
        _WE = WorkflowEnum()

        # initialize logger
        logger = EntryPointLogger()

        # get the input parameters and parse them
        parser = argparse.ArgumentParser(
            description='Implements entry point for the "Icolos" workflow class.'
        )
        parser.add_argument(
            "-conf",
            type=str,
            default=None,
            help="A path to an workflow's configuration file (JSON dictionary) that is to be executed.",
        )
        parser.add_argument(
            "-debug",
            action="store_true",
            help='Set this flag to activate the inbuilt debug logging mode (this will overwrite parameter "-log_conf", if set).',
        )
        parser.add_argument(
            "--global_variables",
            nargs="+",
            default=None,
            type=str,
            help='List of strings, setting global variables with key and value, e.g. "root:/path/to/root".',
        )
        parser.add_argument(
            "--global_settings",
            nargs="+",
            default=None,
            type=str,
            help='List of strings, setting global settings with key and value, e.g. "remove_temporary:False".',
        )
        args, args_unk = parser.parse_known_args()

        if args.conf is None or not os.path.isfile(args.conf):
            raise Exception(
                'Parameter "-conf" must be a relative or absolute path to a configuration (JSON) file.'
            )

        # load configuration
        with open(args.conf) as file:
            conf = file.read().replace("\r", "").replace("\n", "")
            conf = json.loads(conf)

        # set the logging configuration according to parameters
        log_conf = attach_root_path(_LE.PATH_CONFIG_DEFAULT)
        if args.debug:
            log_conf = attach_root_path(_LE.PATH_CONFIG_DEBUG)
        logger = initialize_logging(log_conf_path=log_conf, workflow_conf=conf)

        # update global variables and settings
        conf = parse_header(
            conf=conf,
            args=args,
            entry_point_path=os.path.realpath(__file__),
            logger=logger,
        )

        # generate workflow object
        workflow = WorkFlow(**conf[_WE.WORKFLOW])
        workflow.initialize()

        # execute the whole workflow
        workflow.execute()

        sys.exit(0)


def entry_point():
    IcolosCLI()


if __name__ == "__main__":
    entry_point()
