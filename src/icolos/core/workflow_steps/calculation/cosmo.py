import os
import tempfile
from typing import Tuple, List
from copy import deepcopy

from pydantic import BaseModel

from icolos.utils.execute_external.turbomole import TurbomoleExecutor

from icolos.core.containers.compound import Conformer, Enumeration

from icolos.utils.enums.program_parameters import TurbomoleEnum
from icolos.utils.enums.program_parameters import CosmoOutputEnum
from icolos.utils.enums.compound_enums import ConformerContainerEnum
from icolos.core.workflow_steps.calculation.base import StepCalculationBase
from icolos.core.workflow_steps.step import _LE
from icolos.loggers.logger_utils import log_multiline_string
from icolos.utils.general.files_paths import attach_root_path

_EE = TurbomoleEnum()
_CTE = ConformerContainerEnum()
_COE = CosmoOutputEnum()


class StepCosmo(StepCalculationBase, BaseModel):
    """Step that executes Cosmo.

    Note, that the execution (especially in conjunction with a preceding turbomole step) is relatively complex.
    (1) Take the coord file from the additional data attached to the conformers,
    (2) run Cosmo,
    (3) extract the final XYZ snapshot with x2t,
    (4) translate it to a SDF file with obabel and
    (5) combined the new coordinates with the tags."""

    def __init__(self, **data):
        super().__init__(**data)

        # initialize the executor and test availability
        self._initialize_backend(executor=TurbomoleExecutor)
        self._check_backend_availability()

        # initialize the executor and test availability
        # as they are linked, use a "TurbomoleExecutor" here
        self._initialize_backend(executor=TurbomoleExecutor)
        self._check_backend_availability()

    def _prepare_tmp_input_directories(
        self, enumeration: Enumeration
    ) -> Tuple[List[str], List[str], List[str], List[str]]:
        tmp_dirs = []
        paths_input_cosmofile = []
        paths_config_cosmotherm = []
        paths_output_cosmotherm = []
        for conformer in enumeration:
            # 1) generate all temporary paths
            tmp_dir = tempfile.mkdtemp()
            path_input_cosmofile = os.path.join(tmp_dir, _EE.TM_OUTPUT_COSMOFILE)
            path_config_cosmofile = os.path.join(tmp_dir, _EE.CT_COSMOTHERM_CONFIG_FILE)
            path_output_cosmotherm = os.path.join(
                tmp_dir, _EE.CT_COSMOTHERM_OUTPUT_FILE
            )

            # 2) write-out the COSMO file
            #    Note, that the generation of the COSMO files is part of the Turbomole execution. The reason is, that
            #    the generation is complicated and uses a lot of input form the TM step, thus "cosmoprep" is
            #    executed there.
            if _CTE.EXTRA_DATA_COSMOFILE not in conformer.get_extra_data().keys():
                self._logger.log(
                    f"In order to write out COSMO files, the content needs to be annotated as extra data in the conformers. Have you executed Turbomole before?",
                    _LE.ERROR,
                )
                raise ValueError("Could not find COSMO data to write out - abort.")
            with open(path_input_cosmofile, "w") as f:
                f.writelines(conformer.get_extra_data()[_CTE.EXTRA_DATA_COSMOFILE])

            # 3) add paths
            tmp_dirs.append(tmp_dir)
            paths_input_cosmofile.append(path_input_cosmofile)
            paths_config_cosmotherm.append(path_config_cosmofile)
            paths_output_cosmotherm.append(path_output_cosmotherm)

        return (
            tmp_dirs,
            paths_input_cosmofile,
            paths_config_cosmotherm,
            paths_output_cosmotherm,
        )

    def _execute_run(self, config_path: str):
        result = self._backend_executor.execute(
            command=_EE.CT_COSMOTHERM, arguments=[config_path], check=True
        )
        if _EE.CT_COSMOTHERM_FAIL_STRING in result.stderr:
            self._logger.log(
                f"Execution of {_EE.CT_COSMOTHERM} failed. Error message:", _LE.ERROR
            )
            log_multiline_string(
                logger=self._logger_blank,
                level=_LE.ERROR,
                multi_line_string=result.stdout,
            )

    def _write_config_file(self, config_path: str):
        # by default use the internal configuration, but if one has been specified, use this one
        # note, that the default name of the COSMO file is "mol.cosmo", so this should be used in any config file
        if _EE.CT_CONFIG not in self.settings.arguments.parameters.keys():
            with open(attach_root_path(_EE.CT_CONFIG_DEFAULTPATH), "r") as f:
                config = f.readlines()
            self._logger.log(
                f"Loaded {_EE.CT_COSMOTHERM} configuration from default file {_EE.CT_CONFIG_DEFAULTPATH}.",
                _LE.DEBUG,
            )
        else:
            config = self.settings.arguments.parameters[_EE.CT_CONFIG]
        with open(config_path, "w") as f:
            f.writelines([line.rstrip("\n") + "\n" for line in config])

    def _get_line_by_pattern(self, lines: List[str], pattern: str) -> str:
        for line in lines:
            if pattern in line:
                return line

    def _get_values_from_line(self, line: str) -> List[str]:
        try:
            value_part = line.split(":")[1]
            return value_part.split()
        except Exception:
            return []

    def _annotate_from_output_block(
        self, conformer: Conformer, block: List[str], annotation: dict
    ):
        for key in annotation.keys():
            # get the line with the values
            line = self._get_line_by_pattern(
                lines=block, pattern=annotation[key][_COE.PATTERN]
            )
            if line is None:
                continue

            # get the values and select the one that is to be added
            try:
                values = self._get_values_from_line(line=line)
                value = values[annotation[key][_COE.ELEMENT]]
            except IndexError:
                continue

            # add it as a tag to the conformer; we can replace part of the tag name with e.g. the solvent
            # names if we need to
            conformer.get_molecule().SetProp(key, value)

    def _get_solvents_from_header(self, header: List[str]):
        line_solvents = self._get_line_by_pattern(
            header, pattern=_COE.SOLVENT_BLOCK_HEADER_COMPOUNDS_PATTERN
        )
        return self._get_values_from_line(line_solvents)

    def _get_current_solvent_from_header(self, header: List[str]):
        line_mol_fraction = self._get_line_by_pattern(
            header, pattern=_COE.SOLVENT_BLOCK_HEADER_MOLFRACTION_PATTERN
        )
        solvent_index = self._get_values_from_line(line_mol_fraction).index(
            _COE.SOLVENT_BLOCK_CURRENT_FRACTION_VALUE
        )
        return self._get_solvents_from_header(header)[solvent_index]

    def _parse_general_block(self, lines: List[str], conformer: Conformer):
        general_block = []
        for index in range(len(lines)):
            if _COE.GENERAL_BLOCK_PATTERN_STRING in lines[index]:
                # skip the first lines after the header
                index += 2
                while not lines[index] == "":
                    general_block.append(lines[index])
                    index += 1
                break
        self._annotate_from_output_block(
            conformer=conformer,
            block=general_block,
            annotation=_COE.GENERAL_BLOCK_ANNOTATIONS,
        )

    def _load_solvent_blocks(self, lines: List[str]) -> List[dict]:
        solvent_blocks = []
        index = 0
        while index < len(lines):
            if _COE.SOLVENT_BLOCK_PATTERN_STRING in lines[index]:
                # we need to extract both the header (which solvent?) and the body (actual values)
                new_block = {"header": [], "body": []}
                # go back to start of block
                while (
                    index >= 0 and _COE.SOLVENT_BLOCK_START_PATTERN not in lines[index]
                ):
                    index -= 1

                # extract the header
                while _COE.SOLVENT_BLOCK_BODY_START_PATTERN not in lines[index]:
                    new_block["header"].append(lines[index])
                    index += 1

                # extract the body
                while index < len(lines) and not (
                    lines[index] == "" and lines[index + 1] == ""
                ):
                    new_block["body"].append(lines[index])
                    index += 1
                solvent_blocks.append(new_block)
            index += 1
        return solvent_blocks

    def _annotate_solvent_blocks(
        self, solvent_blocks: List[dict], conformer: Conformer
    ):
        for block_dict in solvent_blocks:
            # get solvent and translate according to internal solvent abbreviation table
            try:
                current_solvent = self._get_current_solvent_from_header(
                    block_dict["header"]
                )
                if current_solvent in _COE.SOLVENT_TRANSLATE_SOLVENT.keys():
                    current_solvent = _COE.SOLVENT_TRANSLATE_SOLVENT[current_solvent]
            except ValueError:
                continue

            # overwrite the solvent name placeholder in and annotate
            template_annotations = deepcopy(_COE.SOLVENT_BLOCK_BODY_ANNOTATIONS)
            annotations = {}
            for key in template_annotations.keys():
                new_key = key.replace(_COE.SOLVENT_REPLACEHOLDER, current_solvent)
                annotations[new_key] = template_annotations[key]

            # annotate
            self._annotate_from_output_block(
                conformer=conformer, block=block_dict["body"], annotation=annotations
            )

    def _parse_output(self, path_output: str, conformer: Conformer):
        # there are two sets of blocks we need to parse: the "general" block, that is always present and, if specified,
        # free energies from solvents ("mixtures")
        # 1) load the file
        with open(path_output, "r") as f:
            lines = f.readlines()
            lines = [line.rstrip("\n") for line in lines]

        # 2) extract the general block: from the match of the pattern line until the second empty line occurs
        #    e.g. "--- Compound 1 (mol) ---\n\nAtomic weights : 111111\n <etc> ...\n\n"
        self._parse_general_block(lines, conformer)

        # 3) extract the solvent blocks (if available)
        #    search for the first occurrence of a Gibb's free energy and expand top until the pattern line is found and
        #    to bottom until more than one empty line is hit; proceed until all blocks are processed
        solvent_blocks = self._load_solvent_blocks(lines)

        if len(solvent_blocks) > 0:
            self._annotate_solvent_blocks(solvent_blocks, conformer)

    def execute(self):
        for compound in self.get_compounds():
            for enumeration in compound.get_enumerations():
                if len(enumeration.get_conformers()) == 0:
                    continue

                # generate copies of the conformers, as to not accidentally manipulate them
                inp_enum = deepcopy(enumeration)

                # prepare the temporary files and retrieve paths (TM config is charge-state dependent!)
                (
                    tmp_dirs,
                    paths_input_cosmofile,
                    paths_config_cosmotherm,
                    paths_output_cosmotherm,
                ) = self._prepare_tmp_input_directories(enumeration=inp_enum)

                # execute individual conformers
                for (
                    tmp_dir,
                    path_config_cosmotherm,
                    conformer,
                    path_output_cosmotherm,
                ) in zip(
                    tmp_dirs,
                    paths_config_cosmotherm,
                    enumeration.get_conformers(),
                    paths_output_cosmotherm,
                ):
                    self._move_to_dir(tmp_dir)

                    # set a necessary environment variable to avoid clashes
                    os.environ[_EE.TM_TURBOTMPDIR] = tmp_dir

                    # write configuration file
                    self._write_config_file(config_path=path_config_cosmotherm)

                    # all ready; start the execution
                    self._execute_run(config_path=path_config_cosmotherm)

                    # parse the results
                    self._parse_output(
                        path_output=path_output_cosmotherm, conformer=conformer
                    )

                # restore working directory and remove temporary files
                self._restore_working_dir()
                for tmp_dir in tmp_dirs:
                    if os.path.isdir(tmp_dir):
                        self._remove_temporary(tmp_dir)

                self._logger.log(
                    f"Executed COSMO for {len(enumeration.get_conformers())} conformers for enumeration {enumeration.get_index_string()}.",
                    _LE.INFO,
                )
