#!/usr/bin/env python
#  coding=utf-8

import os
import sys
import json
import argparse
from datetime import datetime
from icolos.core.composite_agents.workflow import WorkFlow


from icolos.utils.enums.composite_agents_enums import WorkflowEnum
from icolos.utils.enums.logging_enums import LoggingConfigEnum
from icolos.utils.enums.entry_points import ExecutorEnum

from icolos.utils.entry_point_functions.logging_helper_functions import (
    initialize_logging,
)
from icolos.utils.entry_point_functions.parsing_functions import parse_header, log_version_number, get_version_number
from icolos.utils.general.citation_generator import print_citations
from icolos.utils.general.files_paths import attach_root_path


def main():
    # enums
    _LE = LoggingConfigEnum()
    _EE = ExecutorEnum()
    _WE = WorkflowEnum()

    # get the input parameters and parse them
    parser = argparse.ArgumentParser(
        description='Implements entry point for the "Icolos" workflow class.',
        epilog=f"Icolos version: {get_version_number()}"
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
        help="Set this flag to activate the inbuilt debug logging mode.",
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
        help='List of strings, setting global settings with key and value, e.g. "remove_temporary_files:False" "single_directory:True".',
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

    # write the version of the installation used to the logfile
    log_version_number(logger)

    # update global variables and settings
    conf = parse_header(
        conf=conf, args=args, entry_point_path=os.path.realpath(__file__), logger=logger
    )
    # generate workflow object
    workflow = WorkFlow(**conf[_WE.WORKFLOW])
    workflow.initialize()

    # print the logo and ,depending on the backend software used, the citations for the workflow given
    print_citations(workflow._initialized_steps)

    # execute the whole workflow
    st_time = datetime.now()
    workflow.execute()
    exec_time = datetime.now() - st_time
    logger.log(f"Icolos workflow completed. Walltime: {exec_time}.", _LE.INFO)


if __name__ == "__main__":
    main()

    sys.exit(0)
