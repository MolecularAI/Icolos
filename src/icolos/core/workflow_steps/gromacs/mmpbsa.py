import numpy as np
from subprocess import CompletedProcess
import tempfile
from icolos.core.containers.generic import GenericData
from typing import AnyStr, List
from icolos.core.workflow_steps.gromacs.base import StepGromacsBase
from pydantic import BaseModel
from icolos.utils.enums.step_enums import StepGromacsEnum
from icolos.core.workflow_steps.step import _LE
from icolos.utils.enums.program_parameters import GromacsEnum
import os
from icolos.utils.execute_external.gromacs import GromacsExecutor
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

        self._initialize_backend(GromacsExecutor)
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
        threads = self.get_additional_setting(_SGE.THREADS, default=1)
        if threads > 1:
            command = f"mpirun -np {threads} " + command
        self._logger.log(
            f"Executing mmgbsa calculation with {threads} thread(s)", _LE.DEBUG
        )
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

        arguments = ["-f", _SGE.STD_STRUCTURE]
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

    def parse_results(self, wkdirs: List[str]):
        for wkdir in wkdirs:
            results_file = os.path.join(wkdir, "FINAL_RESULTS_MMPBSA.dat")
            with open(results_file, "r") as f:
                lines = f.readlines()
            # for now just parse the final energy value
            result = float([l for l in lines if "DELTA TOTAL" in l][0].split()[2])
            try:
                self.get_topol().properties[_SGE.MMGBSA_DG].append(result)
            except KeyError:
                self.get_topol().properties[_SGE.MMGBSA_DG] = [result]
        avg_dG = np.mean(self.get_topol().properties[_SGE.MMGBSA_DG])
        self._logger.log(f"Average dG = {avg_dG}", _LE.INFO)
        # write the invididual lines to the tmpdir
        parent_dir = os.path.dirname(wkdirs[0])
        with open(os.path.join(parent_dir, "results_summary.dat"), "w") as f:
            for idx, val in enumerate(self.get_topol().properties[_SGE.MMGBSA_DG]):
                f.write(str(idx) + "," + str(val))

    def execute(self) -> None:
        tmp_dir = self._make_tmpdir()
        topol = self.get_topol()
        self._generate_amber_input_file()
        self.write_input_files(tmp_dir, topol=topol)

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
        wkdirs = [
            tempfile.mkdtemp(dir=tmp_dir) for _ in range(len(topol.trajectories.keys()))
        ]

        for i, wkdir in enumerate(wkdirs):
            self._logger.log(
                f"Running MMGBSA on trajectory {i+1} (of {len(wkdirs)})", _LE.DEBUG
            )
            topol.write_structure(wkdir, index=i)
            topol.write_tpr(wkdir, index=i)
            topol.write_trajectory(wkdir, index=i)

            flag_dict = {
                "-i": os.path.join(tmp_dir, _SGE.MMPBSA_IN),
                "-cs": _SGE.STD_TPR,
                "-cg": self._parse_coupling_groups(tmp_dir),
                "-ci": os.path.join(tmp_dir, _SGE.STD_INDEX),
                "-ct": _SGE.STD_XTC,
                "-cp": os.path.join(tmp_dir, _SGE.STD_TOPOL),
                # do not attempt to open the results in the GUI afterwards
                "-nogui": "",
            }

            flag_list = self._parse_arguments(flag_dict=flag_dict)

            result = self._run_mmpbsa(args=flag_list, tmp_dir=wkdir)
        self.parse_results(wkdirs)

        # parse and delete generated output
        self._parse_output(tmp_dir)
        self._remove_temporary(tmp_dir)
