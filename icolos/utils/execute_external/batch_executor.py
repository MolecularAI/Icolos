import os
from shlex import quote
from icolos.utils.execute_external.execute import ExecutorBase
from icolos.utils.enums.program_parameters import SlurmEnum
import subprocess
from typing import List
import time
from tempfile import mkstemp

_SE = SlurmEnum()


class BatchExecutor(ExecutorBase):
    """For execution of batch jobs using either Slurm or SGE scheduler."""

    def __init__(
        self,
        cores: int,
        partition: str,
        time: str,
        mem: str,
        modules: List,
        other_args: dict,
        gres: str,
        prefix_execution=None,
        binary_location=None,
    ):
        super().__init__(
            prefix_execution=prefix_execution, binary_location=binary_location
        )

        self.cores = cores
        self.partition = partition
        self.time = time
        self.mem = mem
        self.modules = modules
        self.other_args = other_args
        self.gres = gres

    def execute(
        self,
        command: str = None,
        arguments: list = None,
        check: bool = True,
        location=None,
        pipe_input=None,
        tmpfile: str = None,
    ):
        if tmpfile is None:
            tmpfile = self.prepare_batch_script(
                command, arguments, pipe_input, location
            )
        sbatch_command = f"sbatch {tmpfile}"
        # execute the batch script
        result = super().execute(
            command=sbatch_command, arguments=[], location=location
        )
        job_id = result.stdout.split()[-1]
        state = self._wait_for_job_completion(job_id=job_id)

        # check the result from slurm
        if check == True:
            if state != _SE.COMPLETED:
                raise subprocess.SubprocessError(
                    f"Subprocess returned non-zero exit status:\n{sbatch_command}\n Status:\n{state}"
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
        raise NotImplementedError(
            "Cannot reliably check, whether a random program executes properly - do not use."
        )

    def _prepare_command(
        self, command: str, arguments: List, pipe_input: str = None
    ) -> str:
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
        return " ".join(complete_command)

    def _wait_for_job_completion(self, job_id):
        completed = False
        state = None
        while completed is False:
            state = self._check_job_status(job_id)
            if state in [_SE.PENDING, _SE.RUNNING]:
                time.sleep(5)
                continue
            elif state == _SE.COMPLETED:
                completed = True
            elif state == _SE.FAILED:
                completed = True

        return state

    def _check_job_status(self, job_id):
        """
        Monitor the status of a previously submitted job, return the result
        """
        command = f"module load slurmtools && jobinfo {job_id}"
        result = subprocess.run(
            command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )

        state = None
        for line in result.stdout.split("\n"):
            if _SE.STATE in line:
                state = line.split(":")[-1].split()[0]
        return state

    def _construct_slurm_header(self):
        header = [
            "#!/bin/bash",
            f"#SBATCH  -c{self.cores}",
            f"#SBATCH -p {self.partition}",
            f"#SBATCH --time={self.time}",
        ]
        header.append(f"#SBATCH --gres={self.gres}")
        for key, value in self.other_args.items():
            header.append(f"#SBATCH {key}={value}")

        for module in self.modules:
            header.append(f"module load {module}")

        return header
