from subprocess import CompletedProcess
from icolos.core.containers.generic import GenericData
from typing import AnyStr, List
from icolos.core.workflow_steps.gromacs.base import StepGromacsBase
from pydantic import BaseModel
from icolos.utils.enums.step_enums import StepGromacsEnum
from icolos.core.workflow_steps.step import _LE
from icolos.utils.enums.program_parameters import GromacsEnum
import os
from icolos.utils.general.files_paths import attach_root_path

_SGE = StepGromacsEnum()
_GE = GromacsEnum()


class StepGMXmmpbsa(StepGromacsBase, BaseModel):
    """
    Execute gmx_MMPBSA, calculates binding free energy of
    protein-ligand complex using single trajectory approximation,
    using Amber's mmpbsa.py script
    """

    def __init__(self, **data):
        super().__init__(**data)

        self._initialize_backend(self._get_gromacs_executor())
        self._check_backend_availability()

    def _get_arg(self, ext) -> AnyStr:
        return self.data.generic.get_argument_by_extension(ext)

    def _generate_amber_input_file(self) -> None:
        input_file = (
            self.settings.additional[_SGE.INPUT_FILE]
            if _SGE.INPUT_FILE in self.settings.additional.keys()
            else None
        )
        # Normally the user should provide an input file to control the mmgbsa protocol
        if input_file is not None and os.path.isfile(input_file):
            self._logger.log(
                f"Using provided AMBER input file at {self.settings.additional[_SGE.INPUT_FILE]}",
                _LE.DEBUG,
            )
            with open(input_file, "r") as f:
                template = GenericData(file_name="mmpbsa.in", file_data=f.read())
        else:
            self._logger.log("No input file found, defaulting to template", _LE.WARNING)
            # parses user arguments and creates the formatted amber input file from the user specification
            with open(attach_root_path(_SGE.DEFAULT_MMPBSA_IN), "r") as f:
                template = GenericData(file_name="mmpbsa.in", file_data=f.read())

        self.data.generic.add_file(template)

    def _parse_arguments(self, flag_dict: dict) -> List:
        args = []
        for flag in self.settings.arguments.flags:
            if flag != "-O":
                args.append(flag)
        for key, value in self.settings.arguments.parameters.items():
            args.append(key)
            args.append(value)
        for key, value in flag_dict.items():
            if key not in args:
                args.append(key)
                args.append(value)

        # capture output
        return args

    def _run_mmpbsa(self, args, tmp_dir) -> CompletedProcess:
        command = _GE.MMPBSA
        self._logger.log(f"Executing mmgbsa calculation in dir {tmp_dir}", _LE.DEBUG)
        result = self._backend_executor.execute(
            command=command, arguments=args, check=True, location=tmp_dir
        )
        for line in result.stdout.split("\n"):
            self._logger_blank.log(line, _LE.INFO)
        for line in result.stderr.split("\n"):
            self._logger_blank.log(line, _LE.INFO)

        return result

    def _parse_coupling_groups(self, tmp_dir) -> AnyStr:
        # parse the coupling groups to their indexes
        output = []
        pipe_input = self.settings.additional[_SGE.COUPLING_GROUPS]

        structure = self.data.generic.get_argument_by_extension(
            _SGE.FIELD_KEY_STRUCTURE
        )
        arguments = ["-f", structure]
        if [f for f in os.listdir(tmp_dir) if f.endswith("ndx")]:
            arguments.extend(["-n", "index.ndx"])
        else:
            arguments.extend(["-o", "index.ndx"])

        result = self._backend_executor.execute(
            command=_GE.MAKE_NDX,
            arguments=arguments,
            location=tmp_dir,
            check=True,
            pipe_input='echo -e "q"',
        )
        for param in pipe_input.split():
            for line in result.stdout.split("\n"):
                parts = line.split()
                if param in line and parts[1] == param:
                    output.append(parts[0])
                    break
        self._logger.log(f"Resolved coupling groups {output}", _LE.DEBUG)
        return " ".join(output)

    def _get_file_from_dir(self, tmp_dir: str, ext: str) -> AnyStr:
        file = [f for f in os.listdir(tmp_dir) if f.endswith(ext)]
        assert len(file) == 1
        return file[0]

    def execute(self) -> None:
        """
        Execute gmx_MMPBSA
        Note: execution using mpirun is not supported for stability reasons
        """
        tmp_dir = self._make_tmpdir()

        self._generate_amber_input_file()
        self._write_input_files(tmp_dir)

        # gmx_MMPBSA requires the coupling groups of the receptor and ligand

        # form any required coupling groups with make_ndx_command before parsing coupling groups
        # e.g. combine protein + cofactor
        ndx_commands = (
            self.settings.additional[_SGE.MAKE_NDX_COMMAND]
            if _SGE.MAKE_NDX_COMMAND in self.settings.additional.keys()
            else None
        )
        if ndx_commands is not None:
            # can run make_ndx multiple times for complex cases, each set of pipe imput must be separated by a semicolon
            for args in ndx_commands.split(";"):
                self._add_index_group(tmp_dir=tmp_dir, pipe_input=args)
        flag_dict = {
            "-i": _SGE.MMPBSA_IN,
            "-cs": self._get_arg("tpr"),
            "-cg": self._parse_coupling_groups(tmp_dir),
            "-ci": self._get_file_from_dir(tmp_dir=tmp_dir, ext="ndx"),
            "-ct": self._get_arg("xtc"),
            "-cp": self._get_arg("top"),
            # do not attempt to open the results in the GUI afterwards
            "-nogui": "",
        }

        flag_list = self._parse_arguments(flag_dict=flag_dict)

        result = self._run_mmpbsa(flag_list, tmp_dir)

        # parse and delete generated output
        self._parse_output(tmp_dir)
        self._remove_temporary(tmp_dir)
