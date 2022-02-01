from copy import deepcopy
from typing import List

from icolos.utils.enums.step_enums import StepBaseEnum, StepFepPlusEnum
from icolos.utils.enums.program_parameters import FepPlusEnum
from icolos.utils.execute_external.fep_plus import FepPlusExecutor

from pydantic import BaseModel, PrivateAttr
import os
from icolos.core.workflow_steps.step import _LE
import time
from icolos.core.workflow_steps.schrodinger.fep_base import StepFEPBase

from icolos.utils.general.icolos_exceptions import StepFailed

_FE = FepPlusEnum()
_SFE = StepFepPlusEnum()
_SBE = StepBaseEnum


class StepFepPlusExec(StepFEPBase, BaseModel):
    """
    Execute the FEP+ workflow, interfaced with AWS
    """

    class Config:
        underscore_attrs_are_private = True

    _job_id = PrivateAttr()

    def __init__(self, **data):
        super().__init__(**data)

        self._initialize_backend(executor=FepPlusExecutor)
        self._check_backend_availability()

        self._job_id = None

    def _parse_arguments(self):
        parameters = deepcopy(self.settings.arguments.parameters)
        arguments = []
        if len(self.settings.arguments.flags) > 0:
            for flag in self.settings.arguments.flags:
                arguments.append(str(flag))
        if parameters:
            for key in parameters.keys():
                arguments.append(key)
                if parameters[key] is not None and parameters[key] != "":
                    arguments.append(str(parameters[key]))
        # for our AWS config, need to set processors per job =1
        if "-ppj" not in arguments:
            arguments.extend(["-ppj", "1"])
            self._logger.log(
                "Set -ppj 1 for AWS execution, since no override was specified",
                _LE.DEBUG,
            )
        if _SFE.RETRIES not in arguments:
            arguments.extend([_SFE.RETRIES, "3"])
        arguments.append(_SFE.FMP_OUTPUT_FILE)

        # remove "-WAIT" if it has been set, as this will interfere with the implementation (and might cause issues
        # due to file system write buffering)
        if _SFE.WAIT_FLAG in arguments:
            self._logger.log(
                "Ignoring -WAIT flag for FEP+ execution (this would interfere with the implementation).",
                _LE.WARNING,
            )
            arguments = [arg for arg in arguments if arg != _SFE.WAIT_FLAG]
        return arguments

    def _unit_test_simulate_output(self, fmp_data, log_data):
        # call this method from the unit instead of the execute method to write out the expected output
        tmp_dir = self._make_tmpdir()
        with open(
            os.path.join(
                tmp_dir,
                f"{self.settings.arguments.parameters[_SFE.JOBNAME_FLAG]}_{_SFE.FMP_OUTPUT_FILE}",
            ),
            "w",
        ) as f:
            f.write(fmp_data)
        with open(
            os.path.join(
                tmp_dir,
                f"{self.settings.arguments.parameters[_SFE.JOBNAME_FLAG]}_{_SFE.LOGFILE}",
            ),
            "w",
        ) as f:
            f.write(fmp_data)
        self._parse_output(tmp_dir)
        self._extract_log_file_data(tmp_dir)
        self._remove_temporary(tmp_dir)

    def _get_job_id(self, result):
        parts = str(result.stdout).split("\n")
        for part in parts:
            if _SFE.JOBID_STRING in part:
                # full_job_id looks something like 549a938d-d2ca-11eb-b9f2-0a6713e9bd3a but only the first part of the
                # hash is needed to access the right job afterwards
                full_job_id = part.split(" ")[1]
                self._job_id = full_job_id.split("-")[0]
                self._logger.log(f"JobId of FEP+ run is {self._job_id}.", _LE.DEBUG)
        if self._job_id is None:
            self._logger.log(
                "Could not obtain JobId after execution - abort.", _LE.ERROR
            )
            raise StepFailed

    def _get_log_file(self) -> List[str]:
        arguments = [
            self._job_id,
            _SFE.FILE_NAME,
            f'{self.settings.arguments.parameters[_SFE.JOBNAME_FLAG]}_{_SFE.LOGFILE}"',
        ]
        logging_result = None
        trials = 0
        while trials < 30000:
            logging_result = self._backend_executor.execute(
                command=_FE.JSC_TAIL_FILE, arguments=arguments, check=False
            )
            if logging_result.returncode == 1:
                time.sleep(30)
                trials += 1
                continue
            elif logging_result.returncode == 0:
                break
        if logging_result is None:
            raise StepFailed("Could not obtain log file from server within time limit.")
        log_lines = str(logging_result.stdout).split("\n")
        return log_lines

    def _get_new_lines(self, old_file) -> List[str]:
        new_lines = self._get_log_file()
        # take the first n lines off the new log file where n is the length of the old log file
        diff = new_lines[len(old_file) - 1 :]
        return diff

    def _wait_for_job_completion(self):
        # get the log file at this state
        log_file = self._get_log_file()
        for line in log_file:
            self._logger_blank.log(line, _LE.INFO)
        # TODO: set maximum (or at least allow to set a maximum)
        while (_SFE.FEP_EXEC_COMPLETE not in log_file) and (
            _SFE.FEP_EXEC_PARTIAL_COMPLETE not in log_file
        ):
            time.sleep(30)
            new_lines = self._get_new_lines(log_file)
            if len(new_lines) > 0:
                for line in new_lines:
                    self._logger_blank.log(line, _LE.INFO)
                    log_file.append(line)

    def _clean_up(self, tmp_dir: str):
        self._remove_temporary(tmp_dir)
        self._job_id = None

    def execute(self):
        # generate the temporary directory and populate it with the required files
        tmp_dir = self._make_tmpdir()
        self.data.generic.write_out_all_files(tmp_dir)

        # check compounds loaded in properly
        if not self.data.compounds:
            self._logger.log(
                f"No compounds were loaded for step {self.step_id}!  If this was intentional you can ignore this warning.",
                _LE.WARNING,
            )

        # obtain the arguments as a list of strings
        arguments = self._parse_arguments()
        self._logger.log(f"Executing FEP+ calculation in {tmp_dir}.", _LE.INFO)

        # execute fep_plus
        self._apply_token_guard()
        result = self._backend_executor.execute(
            command=_FE.FEP_EXECUTOR, arguments=arguments, location=tmp_dir, check=True
        )

        # get job ID from the job server
        self._get_job_id(result)

        # wait for job completion
        self._wait_for_job_completion()

        # extract the edge information from the log file (rather than the annotated map, as this is easier)
        self._parse_output(tmp_dir)
        self._extract_log_file_data(tmp_dir)
        self._logger.log(f"Completed FEP+ execution.", _LE.INFO)

        # clean-up and reset
        self._clean_up(tmp_dir)
