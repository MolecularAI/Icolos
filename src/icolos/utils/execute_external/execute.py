import os
import abc
import subprocess
from shlex import quote

from icolos.utils.enums.logging_enums import LoggingConfigEnum
from icolos.loggers.steplogger import StepLogger

_LE = LoggingConfigEnum()


class ExecutorBase(metaclass=abc.ABCMeta):
    """Virtual base class for the general and program-specific executors."""

    def __init__(self, prefix_execution=None, binary_location=None):
        # if something needs to be attached to the execution string each time, store it here; if not, value is "None"
        self._prefix_execution = prefix_execution
        self._binary_location = binary_location
        # initialise from the step with self.execution.resource dict

    @abc.abstractmethod
    def execute(
        self, command: str, arguments: list, check=True, location=None, pipe_input=None
    ):
        # to avoid security issues, escape the arguments
        arguments = [quote(str(arg)) for arg in arguments]

        # allow for piped input to be passed to binaries
        if pipe_input is not None:
            # pipe_input = self._parse_pipe_input(pipe_input)
            command = pipe_input + " | " + command

        # check, if command (binary) is to be found at a specific location (rather than in $PATH)
        if self._binary_location is not None:
            command = os.path.join(self._binary_location, command)

        # check, if the something needs to be added before the execution of the "rDock" command
        if self._prefix_execution is not None:
            command = self._prefix_execution + " && " + command

        # execute; if "location" is set, change to this directory and execute there
        complete_command = command + " " + " ".join(str(e) for e in arguments)
        complete_command = [complete_command.replace("'", "")]

        old_cwd = os.getcwd()
        if location is not None:
            os.chdir(location)

        # determine whether this is to be run using local resources or as a batch job
        result = subprocess.run(
            complete_command,
            check=False,  # use the manual check to provide better debugginf information than subprocess
            # convert output to string (instead of byte array)
            universal_newlines=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=True,
        )
        # for line in result.stdout.split("\n"):
        #     print(line)
        if check:
            if result.returncode != 0:
                raise subprocess.SubprocessError(
                    f"Subprocess returned non-zero exit status:\n{complete_command}\nReturn code:\n{result.returncode}\nSTDERR:\n{result.stderr}\nSTDOUT:\n{result.stdout}"
                )
        os.chdir(old_cwd)
        return result

    @abc.abstractmethod
    def is_available(self):
        raise NotImplementedError("Overwrite this method in the child class.")


class Executor(ExecutorBase):
    """For execution of command-line programs that do not have any specific executor themselves."""

    def __init__(self, prefix_execution=None, binary_location=None):
        super().__init__(
            prefix_execution=prefix_execution,
            binary_location=binary_location,
        )

    def execute(
        self, command: str, arguments: list, check=True, location=None, pipe_input=None
    ):
        return super().execute(
            command=command,
            arguments=arguments,
            check=check,
            location=location,
            pipe_input=pipe_input,
        )

    def is_available(self):
        raise NotImplementedError(
            "Cannot reliably check, whether a random program executes properly - do not use."
        )


def execution_successful(output: str, success_str: str) -> bool:
    return True if success_str in output else False
