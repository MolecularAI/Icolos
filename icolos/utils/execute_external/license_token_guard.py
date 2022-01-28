import time
from typing import Dict

from pydantic import BaseModel, PrivateAttr

from icolos.utils.execute_external.execute import Executor

from icolos.loggers.steplogger import StepLogger

from icolos.utils.enums.logging_enums import LoggingConfigEnum
from icolos.utils.enums.program_parameters import SchrodingerExecutablesEnum

_EE = SchrodingerExecutablesEnum()
_LE = LoggingConfigEnum()


class TokenGuardParameters(BaseModel):
    prefix_execution: str = None
    binary_location: str = None
    token_pools: Dict
    wait_interval_seconds: int = 30
    wait_limit_seconds: int = 0


class SchrodingerLicenseTokenGuard(BaseModel):
    """Class that checks, whether enough tokens to execute Schrodinger binaries are available."""

    token_guard: TokenGuardParameters

    class Config:
        underscore_attrs_are_private = True

    _logger = PrivateAttr()
    _executor = PrivateAttr()

    def __init__(self, **data):
        super().__init__(**data)

        self._logger = StepLogger()

        # initialize the executor for all "Schrodinger" related calls and also check if it is available
        self._executor = Executor(
            prefix_execution=self.token_guard.prefix_execution,
            binary_location=self.token_guard.binary_location,
        )

    def _get_token_pool_info(self, licadmin_output: list, token_pool: str) -> dict:
        result = {"found": False}
        for line in licadmin_output:
            if token_pool in line:
                parts = line.split(" ")
                if len(parts) == 16:
                    result["total"] = int(parts[6])
                    result["available"] = int(parts[6]) - int(parts[12])
                    result["found"] = True
                break
        return result

    def _check_licstat_output(self, licadmin_output: list) -> bool:
        all_pools_available = True
        for pool_key, pool_token_numbers in self.token_guard.token_pools.items():
            pool_status = self._get_token_pool_info(licadmin_output, pool_key)
            if pool_status["found"]:
                if pool_status["available"] >= pool_token_numbers:
                    self._logger.log(
                        f"Enough tokens available ({pool_status['available']}) to satisfy requirement ({pool_token_numbers} free tokens) for pool {pool_key}.",
                        _LE.DEBUG,
                    )
                else:
                    self._logger.log(
                        f"Not enough tokens available ({pool_status['available']}) to satisfy requirement ({pool_token_numbers} free tokens) for pool {pool_key}.",
                        _LE.DEBUG,
                    )
                    all_pools_available = False
            else:
                all_pools_available = False
                self._logger.log(
                    f"Could not find information on token pool {pool_key}.", _LE.WARNING
                )
        return all_pools_available

    def _get_licstat_output(self):
        result = self._executor.execute(
            command=_EE.LICADMIN, arguments=[_EE.LICADMIN_STAT], check=True
        )
        if result.returncode != 0:
            self._logger.log(
                f"Could not execute the Schrodinger license token guard - do you need to export the licadmin path?",
                _LE.WARNING,
            )
        return result.stdout.split("\n")

    def guard(self) -> bool:
        # loop over the token pools until they are all satisfied or the time limit has run out
        counter = 0
        success = False
        while True:
            if (
                self.token_guard.wait_limit_seconds != 0
                and (counter * self.token_guard.wait_interval_seconds)
                >= self.token_guard.wait_limit_seconds
            ):
                self._logger.log(
                    f"Wait period ({self.token_guard.wait_limit_seconds} seconds) set for Schrodinger token guard has been exceeded.",
                    _LE.ERROR,
                )
                break

            # reload the output from "licadmin"
            # at this stage, the output from licadmin is a list of strings
            licadmin_output = self._get_licstat_output()

            all_pools_available = self._check_licstat_output(
                licadmin_output=licadmin_output
            )
            if all_pools_available:
                self._logger.log(
                    "All token pool requirements for Schrodinger have been met - proceeding.",
                    _LE.DEBUG,
                )
                success = True
                break
            else:
                time.sleep(self.token_guard.wait_interval_seconds)
                counter = counter + 1

        return success
