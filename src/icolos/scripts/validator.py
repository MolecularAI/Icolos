#!/usr/bin/env python
#  coding=utf-8

import os
import sys
import json
import jsonschema
import argparse

from jsonschema import RefResolver

from icolos.config.schemas.schemas import construct_workflow_schema
from icolos.utils.enums.composite_agents_enums import WorkflowEnum
from icolos.utils.enums.logging_enums import LoggingConfigEnum
from icolos.utils.enums.entry_points import ExecutorEnum

from icolos.utils.entry_point_functions.parsing_functions import (
    get_version_number,
    version_match,
    get_versions_as_strings,
    get_config_version_number,
)
from icolos.utils.general.citation_generator import ConsoleColours

# enums
_LE = LoggingConfigEnum()
_EE = ExecutorEnum()
_WE = WorkflowEnum()

_CC = ConsoleColours()


def main():
    reports = {_LE.INFO: 0, _LE.WARNING: 0, _LE.DEBUG: 0, _LE.ERROR: 0}

    def _report(message: str, level: str = _LE.INFO):
        levels_prefix = {
            _LE.INFO: f"{ConsoleColours.GREEN}INFO:{ConsoleColours.ENDC}",
            _LE.WARNING: f"{ConsoleColours.ORANGE}WARNING:{ConsoleColours.ENDC}",
            _LE.ERROR: f"{ConsoleColours.FAIL}{ConsoleColours.BOLD}ERROR:{ConsoleColours.ENDC}",
            _LE.DEBUG: f"{ConsoleColours.HEADER}DEBUG:{ConsoleColours.ENDC}",
        }
        reports[level] += 1
        print(levels_prefix[level], message, "\n")

    def _summary():
        print("------------------------\nValidation summary\n------------------------")
        print(f"# errors:", reports[_LE.ERROR])
        print(f"# warnings:", reports[_LE.WARNING])
        print("------------------------")

        # if an error occurred, return a failure code (otherwise, even if warnings happened, return 0)
        if reports[_LE.ERROR] > 0:
            print(
                f"\nJSON status: {ConsoleColours.FAIL}{ConsoleColours.BOLD}INVALID{ConsoleColours.ENDC}"
            )
            return 1
        else:
            print(
                f"\nJSON status: {ConsoleColours.GREEN}{ConsoleColours.BOLD}VALID{ConsoleColours.ENDC}"
            )
            return 0

    # get the input parameters and parse them
    version = get_version_number()
    parser = argparse.ArgumentParser(
        description='Implements entry point for the "Icolos validator", checking a JSON configuration file for issues.',
        epilog=f"Icolos version: {'unknown' if version is None else version}",
    )
    parser.add_argument(
        "-conf",
        type=str,
        default=None,
        help="A path to an workflow's configuration file (JSON dictionary) that is to be validated.",
    )
    args, args_unk = parser.parse_known_args()

    if args.conf is None or not os.path.isfile(args.conf):
        raise Exception(
            'Parameter "-conf" must be a relative or absolute path to a configuration (JSON) file.'
        )

    # load configuration
    with open(args.conf) as file:
        conf = file.read().replace("\r", "").replace("\n", "")
        try:
            conf = json.loads(conf)
        except json.decoder.JSONDecodeError as e:
            _report(
                f"Formatting error (missing quotation marks, invalid values, ...) in JSON, error message:\n{str(e)}. Aborting.",
                _LE.ERROR,
            )
            return _summary()

    # check version matching
    version_config = get_config_version_number(conf)
    if version_config is None:
        _report(
            f'Version number of configuration not specified, add to "header" block.',
            _LE.WARNING,
        )
    if not version_match(conf):
        version_config, version_installation = get_versions_as_strings(conf)
        _report(
            f"Versions of installation ({version_installation}) and configuration ({version_config}) do not match or are not specified.",
            _LE.WARNING,
        )

    # load the sub-schemas (e.g. for the header region and the steps etc.), construct a schema for the entire workflow
    # and validate the input JSON against it
    workflow_schema, schema_dir = construct_workflow_schema()
    schema_path = "file:///{0}/".format(schema_dir.replace("\\", "/"))
    resolver = RefResolver(schema_path, workflow_schema)
    validator = jsonschema.Draft7Validator(workflow_schema, resolver)
    for e in validator.iter_errors(conf):
        _report(str(e), level=_LE.ERROR)

    # TODO: add checks on specified paths (if not matched, give warning)
    # TODO: check global variables and settings (if not provided, give warning)

    return _summary()


if __name__ == "__main__":
    ret = main()
    sys.exit(ret)
