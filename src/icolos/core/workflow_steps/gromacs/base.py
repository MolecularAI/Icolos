from icolos.core.containers.generic import GenericData
from icolos.core.containers.gmx_state import GromacsState
from icolos.utils.enums.step_enums import StepGromacsEnum
from pydantic import BaseModel
import os
from typing import List
from icolos.core.workflow_steps.step import StepBase
from icolos.core.workflow_steps.step import _LE
import re
from copy import deepcopy
from distutils.dir_util import copy_tree
from icolos.utils.enums.program_parameters import GromacsEnum
import dill

_SGE = StepGromacsEnum()
_GE = GromacsEnum()


class StepGromacsBase(StepBase, BaseModel):
    def __init__(self, **data):
        super().__init__(**data)

    def write_input_files(self, tmp_dir: str, topol: GromacsState = None):
        """
        Defaults to writing all available data from the workflow's topology object, unless that file has been specified in generic data, which will be prioritised
        """

        # Normally this should be handled by setting GMXLIB env variable, but for some programs (gmx_MMPBSA), this doesn't work and non-standard forcefields
        # need to be in the working directory
        if _SGE.FORCEFIELD in self.settings.additional:
            copy_tree(
                self.settings.additional[_SGE.FORCEFIELD],
                os.path.join(
                    tmp_dir, self.settings.additional[_SGE.FORCEFIELD].split("/")[-1]
                ),
            )
            self._logger.log(
                f"Copied forcefield at {self.settings.additional[_SGE.FORCEFIELD]} to the working "
                f"directory at {tmp_dir}",
                _LE.INFO,
            )

        self._logger.log(
            f"Writing input files to working directory at {tmp_dir}", _LE.DEBUG
        )
        # if we have specified a topol object, write out these files unless overwritten by an argument from the config
        if topol is not None:
            # handle structure files
            if not self.data.generic.get_files_by_extension("gro"):
                if topol.structures:
                    topol.write_structure(tmp_dir)
            if not self.data.generic.get_files_by_extension("top"):
                if topol.top_lines:
                    topol.write_topol(tmp_dir)
            if not self.data.generic.get_files_by_extension("tpr"):
                if topol.tprs:
                    topol.write_tpr(tmp_dir)
            if not self.data.generic.get_files_by_extension("xtc"):
                if topol.trajectories:
                    topol.write_trajectory(tmp_dir)
            if not self.data.generic.get_files_by_extension("ndx"):
                if topol.ndx:
                    topol.write_ndx(tmp_dir)

        for file in self.data.generic.get_flattened_files():
            file.write(tmp_dir)

    def pickle_topol(self, topol: GromacsState, tmp_dir: str):
        """
        Write the Icolos internal state to a pickled file in the tmpdir
        """
        with open(os.path.join(tmp_dir, "icolos_topol_state.pkl"), "wb") as f:
            dill.dump(topol, f)

    def load_topol(self, file: GenericData) -> GromacsState:
        data = file.get_data()
        return dill.load(data)

    def _parse_arguments(self, flag_dict: dict, args: list = None) -> List:
        arguments = args if args is not None else []
        # first add the settings from the command line
        for key in self.settings.arguments.parameters.keys():
            arguments.append(key)
            arguments.append(str(self.settings.arguments.parameters[key]))
        for flag in self.settings.arguments.flags:
            arguments.append(str(flag))
        for key, value in flag_dict.items():
            # only add defaults if they have not been specified in the json
            if key not in arguments:
                arguments.append(key)
                arguments.append(value)
        return arguments

    def _copy_fields_dict(self):
        try:
            update_dictionary = deepcopy(self.settings.additional[_SGE.FIELDS])
            return update_dictionary
        except KeyError:
            self._logger.log(
                "Update dictionary not present, will use provided mdp file without further modification",
                _LE.WARNING,
            )
            return {}

    def generate_output_file(self, in_file):
        parts = in_file.split(".")
        return parts[0] + "_out" + "." + parts[1]

    def _modify_config_file(
        self, tmp_dir: str, config_file: GenericData, update_dict: dict
    ):
        file_data = config_file.get_data()
        for key, value in update_dict.items():
            pattern = rf"({key})(\s*=\s*)[a-zA-Z0-9\s\_]*(\s*;)"
            pattern = re.compile(pattern)
            matches = re.findall(pattern, file_data)
            if len(matches) == 0:
                self._logger.log(
                    f"Specified key {key} was not found in the mdp file, value was not changed!",
                    _LE.WARNING,
                )
            else:

                file_data = re.sub(pattern, rf"\1\2 {value} \3", file_data)
                self._logger.log(
                    f"Replaced field {key} of mdp file with value {value}", _LE.DEBUG
                )
        self._logger.log(f"Final MDP file for step {self.step_id} is: ", _LE.DEBUG)
        for line in file_data.split("\n"):
            self._logger_blank.log(line, _LE.DEBUG)
        config_file.set_data(file_data)
        config_file.write(tmp_dir)

    def _generate_index_groups(self, tmp_dir: str):
        # dump the first structure file to the tmpdir
        structure = self.get_topol().structures[0]
        structure.write(tmp_dir)
        args = ["-f", structure.get_file_name()]
        ndx_list = [f for f in os.listdir(tmp_dir) if f.endswith(_SGE.FIELD_KEY_NDX)]
        if len(ndx_list) == 1:
            args.extend(["-n", ndx_list[0]])
        result = self._backend_executor.execute(
            command=_GE.MAKE_NDX,
            arguments=args,
            location=tmp_dir,
            check=True,
            pipe_input='echo -e "q"',
        )
        self.get_topol().set_ndx(tmp_dir)
        return result

    def construct_pipe_arguments(self, tmp_dir, params) -> str:
        """
        Constructs the pipe arguments to be passed to gromacs interactive programs
        """
        # look up the groups that have been passed, try to identify the group number in the corresponding index file

        result = self._generate_index_groups(tmp_dir)
        output = ['echo -e "']
        for param in params.split():
            if param == "or":
                output.append("|")
            elif param == "and":
                output.append("&")
            elif param == "not":
                output.append("!")
            elif param == ";":
                output.append("\n")
            else:
                added_one = False
                for line in result.stdout.split("\n"):
                    parts = line.split()
                    if param in parts and param == parts[1]:
                        idx = parts[0]
                        # print("found index", idx, f"for {param}")
                        added_one = True
                        output.append(idx)
                        break
                if not added_one:
                    output.append(param)
        output.append('\nq"')
        self._logger.log(f"Constructed pipe input {' '.join(output)}", _LE.DEBUG)
        return " ".join(output)

    def _add_index_group(self, tmp_dir, pipe_input):
        ndx_args = [
            "-f",
            os.path.join(tmp_dir, _SGE.STD_STRUCTURE),
            "-o",
            os.path.join(tmp_dir, _SGE.STD_INDEX),
        ]
        self._logger.log(
            f"Added group to index file using command {pipe_input}",
            _LE.DEBUG,
        )
        result = self._backend_executor.execute(
            command=_GE.MAKE_NDX,
            arguments=ndx_args,
            location=tmp_dir,
            check=True,
            pipe_input=self.construct_pipe_arguments(tmp_dir, pipe_input),
        )
        for line in result.stdout.split("\n"):
            self._logger_blank.log(line, _LE.INFO)
        self.get_topol().set_ndx(tmp_dir)
