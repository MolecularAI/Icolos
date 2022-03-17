import os
from shlex import quote
from icolos.loggers.steplogger import StepLogger
from icolos.utils.enums.logging_enums import LoggingConfigEnum
from icolos.utils.execute_external.execute import ExecutorBase
from icolos.utils.enums.program_parameters import SlurmEnum
import subprocess
from typing import List
import time
from tempfile import mkstemp

_SE = SlurmEnum()
logger = StepLogger()
_LE = LoggingConfigEnum()


class SlurmExecutor(ExecutorBase):
    """For execution of batch jobs to a Slurm cluster."""

    def __init__(
        self,
        cores: int,
        partition: str,
        time: str,
        mem: str,
        tasks: str,
        modules: List,
        other_args: dict,
        gres: str,
        additional_lines: List,
        prefix_execution=None,
        binary_location=None,
    ):
        super().__init__(prefix_execution=None, binary_location=None)

        self.cores = cores
        self.partition = partition
        self.time = time
        self.mem = mem
        self.tasks = tasks
        self.modules = modules
        self.other_args = other_args
        self.gres = gres
        self.additional_lines = additional_lines
        self._script_prefix_execution = prefix_execution
        self._script_binary_location = binary_location

    def execute(
        self,
        command: str = None,
        arguments: list = None,
        check: bool = True,
        location=None,
        pipe_input=None,
        tmpfile: str = None,
    ):
        """
        Creates and executes the batch script using the provided resource requirements
        If a path to an existing batch script has not been passed via tmpfile,  it is created
        Attempts to sbatch the jobscript, falling back on bash to provide compatibility with workstations (in this case #SLURM lines are ignored, and the execution becomes blocking)
        """
        if tmpfile is None:
            tmpfile = self.prepare_batch_script(
                command, arguments, pipe_input, location
            )
        if self.is_available():
            launch_command = f"sbatch {tmpfile}"
        else:
            logger.log(
                "Warning - Slurm was not found, falling back to local execution!",
                _LE.WARNING,
            )
            launch_command = f"bash {tmpfile}"
        # execute the batch script
        result = super().execute(
            command=launch_command, arguments=[], location=location, check=check
        )

        # either monitor the job id, or resort to parsing the log file
        if self.is_available():
            job_id = result.stdout.split()[-1]
            state = self._wait_for_job_completion(job_id=job_id)
        # if using local resources, bash call is blocking, no need to monitor, just wait for result to return
        else:
            state = _SE.COMPLETED if result.returncode == 0 else _SE.FAILED

        # check the result from slurm
        if check == True:
            if state != _SE.COMPLETED:
                raise subprocess.SubprocessError(
                    f"Subprocess returned non-zero exit status:\n{launch_command}\n Status:\n{state}"
                )
        return state

    def prepare_batch_script(
        self, command, arguments: List, pipe_input: str = None, location=None
    ):
        """
        Write a batch script to the specified location
        """
        batch_script = self._construct_slurm_header()
        command = self._prepare_command(command, arguments, pipe_input)
        if isinstance(command, str):
            command = [command]
        for cmd in command:
            batch_script.append(cmd)
        _, tmpfile = mkstemp(dir=location, suffix=".sh")
        with open(tmpfile, "w") as f:
            for line in batch_script:
                f.write(line)
                f.write("\n")

        return tmpfile

    def is_available(self):
        command = "sbatch -h"
        result = super().execute(command=command, arguments=[], check=False)
        if any(["Usage: sbatch" in l for l in result.stdout.split("\n")]):
            return True
        return False

    def _prepare_command(
        self, command: str, arguments: List, pipe_input: str = None
    ) -> str:
        arguments = [quote(str(arg)) for arg in arguments]

        # allow for piped input to be passed to binaries
        if pipe_input is not None:
            # pipe_input = self._parse_pipe_input(pipe_input)
            command = pipe_input + " | " + command

        # check, if command (binary) is to be found at a specific location (rather than in $PATH)
        if self._script_binary_location is not None:
            command = os.path.join(self._script_binary_location, command)

        # check, if the something needs to be added before the execution of the "rDock" command
        if self._prefix_execution is not None:
            command = self._script_prefix_execution + " && " + command

        # execute; if "location" is set, change to this directory and execute there
        complete_command = command + " " + " ".join(str(e) for e in arguments)
        complete_command = [complete_command.replace("'", "")]
        return " ".join(complete_command)

    def _wait_for_job_completion(self, job_id: str):
        completed = False
        state = None
        logger.log(f"Monitoring slurm job {job_id}", _LE.DEBUG)
        while completed is False:
            state = self._check_job_status(job_id)
            if state in [_SE.PENDING, _SE.RUNNING, None]:
                time.sleep(60)
                continue
            elif state == _SE.COMPLETED:
                completed = True
            elif state == _SE.FAILED:
                completed = True
        return state

    def _tail_log_file(
        self,
        location: str,
        completed_line: str = "Finished mdrun",
        failed_line: str = "Fatal error",
    ):
        # TODO, this is not very robust at the moment!
        completed = False
        state = None
        while not completed:
            try:
                with open(os.path.join(location, "md.log"), "r") as f:
                    lines = f.readlines()
                completed = any([completed_line in l for l in lines])
                state = _SE.COMPLETED
                failed = any([failed_line in l for l in lines])
                if failed:
                    state = _SE.FAILED
                    for line in lines[-40:]:
                        print(line)
                    completed = True
            except FileNotFoundError:
                logger.log("log file not found, sleeping", _LE.DEBUG)
                time.sleep(10)
        return state

    def _check_job_status(self, job_id):
        """
        Monitor the status of a previously submitted job, return the result
        """
        # use the entrypoint included in the Icolos install

        command = f"sacct -j {job_id} --parsable --noheader -a"
        result = subprocess.run(
            command,
            shell=True,
            universal_newlines=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if result.stdout:
            state = result.stdout.split("\n")[0].split("|")[5]
        else:
            state = None
        return state

    def _construct_slurm_header(self):
        header = [
            "#!/bin/bash",
            f"#SBATCH --partition={self.partition}",
            f"#SBATCH --time={self.time}",
        ]
        if self.mem is not None:
            header.append(f"#SBATCH --mem={self.mem}")
        if self.cores is not None:
            header.append(f"#SBATCH -c{self.cores}")
        if self.tasks is not None:
            header.append(f"#SBATCH --tasks={self.tasks}")
        if self.gres is not None:
            header.append(f"#SBATCH --gres={self.gres}")
        for key, value in self.other_args.items():
            header.append(f"#SBATCH {key}={value}")

        for module in self.modules:
            header.append(f"module load {module}")

        # add any other specified lines
        for line in self.additional_lines:
            header.append(line)

        return header
