import os
import glob
import tempfile
import time
from typing import List, Tuple

from pydantic import BaseModel, Field
from rdkit import Chem
from copy import deepcopy

from icolos.utils.enums.step_enums import StepBaseEnum, StepGoldEnum
from icolos.utils.execute_external.gold import GoldExecutor
from icolos.utils.general.files_paths import gen_tmp_file
from icolos.core.containers.compound import Conformer
from icolos.utils.enums.program_parameters import GoldOutputEnum, GoldEnum
from icolos.core.workflow_steps.step import _LE, StepBase
from icolos.utils.general.parallelization import Subtask, SubtaskContainer, Parallelizer
from icolos.core.step_utils.gold_config import (
    ConfigAutomaticSettings,
    ConfigPopulation,
    ConfigGeneticOperators,
    ConfigFloodFill,
    ConfigDataFiles,
    ConfigFlags,
    ConfigProteinData,
    ConfigFitnessFunctionSettings,
    ConfigSaveOptions,
    ConfigCovalentBonding,
    ConfigConstraints,
    ConfigTermination,
)

_SBE = StepBaseEnum
_CGE = GoldEnum()
_SGE = StepGoldEnum()
_GOE = GoldOutputEnum()


class GoldConfiguration(BaseModel):
    automatic_settings: ConfigAutomaticSettings = Field(
        alias="AUTOMATIC SETTINGS", default=ConfigAutomaticSettings()
    )
    population: ConfigPopulation = Field(alias="POPULATION", default=ConfigPopulation())
    genetic_operators: ConfigGeneticOperators = Field(
        alias="GENETIC OPERATORS", default=ConfigGeneticOperators()
    )
    flood_fill: ConfigFloodFill = Field(alias="FLOOD FILL", default=None)
    data_files: ConfigDataFiles = Field(alias="DATA FILES", default=ConfigDataFiles())
    flags: ConfigFlags = Field(alias="FLAGS", default=ConfigFlags())
    termination: ConfigTermination = Field(
        alias="TERMINATION", default=ConfigTermination()
    )
    constraints: ConfigConstraints = Field(
        alias="CONSTRAINTS", default=ConfigConstraints()
    )
    covalent_bonding: ConfigCovalentBonding = Field(
        alias="COVALENT BONDING", default=ConfigCovalentBonding()
    )
    save_options: ConfigSaveOptions = Field(
        alias="SAVE OPTIONS", default=ConfigSaveOptions()
    )
    fitness_function_settings: ConfigFitnessFunctionSettings = Field(
        alias="FITNESS FUNCTION SETTINGS", default=ConfigFitnessFunctionSettings()
    )
    protein_data: ConfigProteinData = Field(alias="PROTEIN DATA", default=None)


class GoldAdditional(BaseModel):
    configuration: GoldConfiguration = GoldConfiguration()
    gold_config_file: str = None
    docking_score_tag: str = "Gold.Goldscore.Fitness"
    grid_ids: List[str] = ["grid0"]


class StepGold(StepBase, BaseModel):

    gold_additional: GoldAdditional = None
    _gold_executor: GoldExecutor = None

    def __init__(self, **data):
        super().__init__(**data)

        # initialize the executor and test availability
        self._initialize_backend(executor=GoldExecutor)
        self._check_backend_availability()

        # set Gold specific settings
        self.gold_additional = GoldAdditional(**self.settings.additional)

    def _set_docking_score(self, conformer: Chem.Mol) -> bool:
        tag = self.gold_additional.docking_score_tag
        try:
            docking_score = conformer.GetProp(tag)
        except KeyError:
            self._logger.log(
                f'Could not find tag "{tag}" in conformer, use "docking_score_tag" to adjust.',
                _LE.DEBUG,
            )
            return False
        conformer.SetProp(_SBE.ANNOTATION_TAG_DOCKING_SCORE, str(docking_score))
        return True

    def _generate_temporary_input_output_files(
        self, batch: List[List[Subtask]]
    ) -> Tuple[List[str], List[List[str]], List[List[str]]]:
        tmp_output_dirs = []
        tmp_input_sdf_paths = []
        tmp_output_sdf_paths = []

        for next_subtask_list in batch:
            # generate temporary input files and output directory
            cur_tmp_output_dir = tempfile.mkdtemp()

            # write-out the temporary input files: one molecule per file each
            one_written = False
            list_input_files = []
            list_output_files = []
            for subtask in next_subtask_list:
                enumeration = subtask.data
                mol = deepcopy(enumeration.get_molecule())
                if mol is not None:
                    _, cur_tmp_sdf = gen_tmp_file(suffix=".sdf", dir=cur_tmp_output_dir)
                    writer = Chem.SDWriter(cur_tmp_sdf)
                    mol.SetProp("_Name", enumeration.get_index_string())
                    writer.write(mol)
                    writer.close()
                    list_input_files.append(cur_tmp_sdf)

                    # output files for gold look something like: ranked_tmp1cwtavoh_m1_2.sdf
                    # where: (1) "ranked_" is a constant, (2) followed by the input file name,
                    #        (3) followed by the molecule # (here always "m1" as only 1 compound/file at the moment),
                    #        (4) followed by the pose # and finally (5) the ".sdf" ending
                    # use '*' for getting all files
                    basename_filepart = os.path.basename(cur_tmp_sdf).rstrip(".sdf")
                    list_output_files.append(
                        os.path.join(
                            cur_tmp_output_dir,
                            "".join(["ranked_", basename_filepart, "_m1_*.sdf"]),
                        )
                    )
                    one_written = True
            if one_written is False:
                # no compound left from this batch: remove the temporary folder and skip this one
                self._remove_temporary(cur_tmp_output_dir)
                continue

            # since all went well, add input and output paths here
            tmp_input_sdf_paths.append(list_input_files)
            tmp_output_sdf_paths.append(list_output_files)
            tmp_output_dirs.append(cur_tmp_output_dir)
        return (tmp_output_dirs, tmp_input_sdf_paths, tmp_output_sdf_paths)

    def _execute_gold(self):
        # get number of sublists in batch and initialize Parallelizer
        gold_parallelizer = Parallelizer(func=self._run_subjob)

        # continue until everything is successfully done or number of retries have been exceeded
        while self._subtask_container.done() is False:
            next_batch = self._get_sublists(get_first_n_lists=self._get_number_cores())

            # generate paths and initialize molecules (so that if they fail, this can be covered)
            (
                tmp_output_dirs,
                tmp_input_sdf_paths,
                tmp_output_sdf_paths,
            ) = self._generate_temporary_input_output_files(next_batch)

            # execute the current batch in parallel; hand over lists of parameters (will be handled by Parallelizer)
            # also increment the tries and set the status to "failed" (don't do that inside subprocess, as data is
            # copied, not shared!)
            _ = [sub.increment_tries() for element in next_batch for sub in element]
            _ = [sub.set_status_failed() for element in next_batch for sub in element]
            gold_parallelizer.execute_parallel(
                output_dir=tmp_output_dirs, input_paths=tmp_input_sdf_paths
            )

            # parse the output of that particular batch and remove temporary files
            self._parse_gold_output(
                tmp_output_paths=tmp_output_sdf_paths, batch=next_batch
            )

            # clean-up
            self._remove_temporary(tmp_output_dirs)

            # print the progress for this execution
            self._log_execution_progress()

    def _parse_gold_output(
        self, tmp_output_paths: List[List[str]], batch: List[List[Subtask]]
    ):
        def _update_subtask(sublist: List[Subtask], enum_identifier: str):
            for task in sublist:
                if task.data.get_index_string() == enum_identifier:
                    task.set_status_success()

        grid_id = self.gold_additional.grid_ids[0]
        grid_path = self.gold_additional.configuration.protein_data.protein_datafile

        for batch_id in range(len(batch)):
            cur_sublist = batch[batch_id]
            list_comp_output = tmp_output_paths[batch_id]

            for compound_number in range(len(list_comp_output)):

                # expand the output paths to cover all pose files for this compound
                # one input sdf matches to x output sdfs (for x poses)
                comp_sdf_output_paths = glob.glob(list_comp_output[compound_number])

                # loop over all output files
                for comp_output in comp_sdf_output_paths:
                    # this is a protection against the case where empty (file size == 0 bytes) files
                    # are generated due to a failure during docking
                    if (
                        not os.path.isfile(comp_output)
                        or os.path.getsize(comp_output) == 0
                    ):
                        continue

                    mol_supplier = Chem.SDMolSupplier(comp_output, removeHs=False)
                    for mol in mol_supplier:
                        if mol is None:
                            continue

                        # The name is something like: 10:0|tmp8x3rtg4d|sdf|1|dock7
                        cur_enumeration_name = str(mol.GetProp("_Name")).split("|")[0]

                        # add the information on the actual grid used
                        mol.SetProp(_SBE.ANNOTATION_GRID_ID, str(grid_id))
                        mol.SetProp(_SBE.ANNOTATION_GRID_PATH, str(grid_path))
                        mol.SetProp(
                            _SBE.ANNOTATION_GRID_FILENAME, os.path.basename(grid_path)
                        )

                        # if no docking score is attached (i.e. the molecule is a receptor or so, skip it)
                        if self._set_docking_score(mol) is not True:
                            continue

                        # add molecule to the appropriate ligand
                        for compound in self.get_compounds():
                            for enumeration in compound:
                                if (
                                    enumeration.get_index_string()
                                    == cur_enumeration_name
                                ):
                                    new_conformer = Conformer(
                                        conformer=mol,
                                        conformer_id=None,
                                        enumeration_object=enumeration,
                                    )
                                    enumeration.add_conformer(
                                        new_conformer, auto_update=True
                                    )
                                    _update_subtask(
                                        cur_sublist,
                                        enum_identifier=cur_enumeration_name,
                                    )
                                    break

    def _run_subjob(self, output_dir: str, input_paths: List[str]):
        # generate the appropriate config; note that input and output paths are referring to one compound each
        config_path = os.path.join(output_dir, "gold.config")
        self.generate_config_file(path=config_path, ligand_files=input_paths)

        # set up arguments list and execute; change path to temporary sub-folder to avoid cluttering with files
        arguments = [config_path]
        old_dir = os.getcwd()
        os.chdir(output_dir)
        execution_result = self._backend_executor.execute(
            command=_CGE.GOLD_CALL, arguments=arguments, check=True
        )
        os.chdir(old_dir)
        time.sleep(3)

    def _config_from_file(self, ligand_files: List[str]) -> List[str]:
        config = []
        with open(self.gold_additional.gold_config_file, "r") as f:
            buffer = f.readlines()
            idx = 0
            while idx < len(buffer):
                line = buffer[idx].rstrip("\n")
                config.append(line)
                if _SGE.DATA_FILES in line:
                    # skip over all defined ligand files
                    while _SGE.LIGAND_DATA_FILE in buffer[idx + 1]:
                        self._logger.log(
                            'Skipping pre-defined ligand in configuration file (remove lines starting with "ligand_data_file" to suppress this information).',
                            _LE.DEBUG,
                        )
                        idx = idx + 1

                    # add ligand files as required
                    for ligand_file in ligand_files:
                        config.append(
                            " ".join([_SGE.LIGAND_DATA_FILE, ligand_file, "10"])
                        )
                idx = idx + 1
        return config

    def _config_from_default(self, ligand_files: List[str]) -> List[str]:
        def block_indent(block_name: str) -> str:
            return "".join([_SGE.BLOCK_INDENT, block_name])

        def empty_line() -> str:
            return ""

        def add_block(list_lines: List[str], block: dict):
            for key, value in block.items():
                if key == _SGE.LIGAND_DATA_FILE:
                    # this needs to be overwritten by the actual compound collection
                    self._logger.log(
                        f"Do not set ligand data files explicitly, use the compound handover - skipping.",
                        _LE.WARNING,
                    )
                    continue
                list_lines.append(" ".join([key, "=", str(value)]))
            list_lines.append(empty_line())

        config = [block_indent(_SGE.CONFIGURATION_START), empty_line()]

        # automatic settings
        config.append(block_indent(_SGE.AUTOMATIC_SETTINGS))
        add_block(config, self.gold_additional.configuration.automatic_settings.dict())

        # population
        config.append(block_indent(_SGE.POPULATION))
        add_block(config, self.gold_additional.configuration.population.dict())

        # genetic operators
        config.append(block_indent(_SGE.GENETIC_OPERATORS))
        add_block(config, self.gold_additional.configuration.genetic_operators.dict())

        # flood fill
        config.append(block_indent(_SGE.FLOOD_FILL))
        add_block(config, self.gold_additional.configuration.flood_fill.dict())

        # data files
        config.append(block_indent(_SGE.DATA_FILES))
        for ligand_file in ligand_files:
            config.append(" ".join([_SGE.LIGAND_DATA_FILE, ligand_file, "10"]))
        add_block(config, self.gold_additional.configuration.data_files.dict())

        # flags
        config.append(block_indent(_SGE.FLAGS))
        add_block(config, self.gold_additional.configuration.flags.dict())

        # termination
        config.append(block_indent(_SGE.TERMINATION))
        add_block(config, self.gold_additional.configuration.termination.dict())

        # constraints
        config.append(block_indent(_SGE.CONSTRAINTS))
        add_block(config, self.gold_additional.configuration.constraints.dict())

        # covalent bonding
        config.append(block_indent(_SGE.COVALENT_BONDING))
        add_block(config, self.gold_additional.configuration.covalent_bonding.dict())

        # save options
        config.append(block_indent(_SGE.SAVE_OPTIONS))
        add_block(config, self.gold_additional.configuration.save_options.dict())

        # fitness function settings
        config.append(block_indent(_SGE.FITNESS_FUNCTION_SETTINGS))
        add_block(
            config, self.gold_additional.configuration.fitness_function_settings.dict()
        )

        # protein data
        config.append(block_indent(_SGE.PROTEIN_DATA))
        add_block(config, self.gold_additional.configuration.protein_data.dict())

        return config

    def generate_config_file(self, path: str, ligand_files: List[str]):
        # very complicated custom format for GOLD; see file in "external_documentation/gold.conf" for details

        if self.gold_additional.gold_config_file is None:
            # generate config file step by step from default and overwrite set elements
            # TODO: constraints are not yet supported (and probably a few other fine-grained settings)
            config = self._config_from_default(ligand_files=ligand_files)
        else:
            # load path to config file and replace "ligand_files" as required; as configuration element is ignored in
            # this mode, check that it is not set
            if _SGE.CONFIGURATION in self.gold_additional.dict().keys():
                self._logger.log(
                    'The "configuration" elements are ignored when a config file path is specified.',
                    _LE.WARNING,
                )
            config = self._config_from_file(ligand_files=ligand_files)

        # write out
        with open(path, "w") as f:
            for line in config:
                f.write(line + "\n")

    def execute(self):
        # Note: This step only supports one grid at a time, ensemble docking is taken care of at the workflow level!

        # in order to be able to efficiently execute Gold on the enumeration level, all of them have to be unrolled
        # Note: As they retain their respective Compound object, the attribution later on is simple
        all_enumerations = []
        for compound in self.get_compounds():
            all_enumerations = all_enumerations + compound.get_enumerations()
            for enumeration in compound:
                enumeration.clear_conformers()

        # split into sublists, according to the settings
        self._subtask_container = SubtaskContainer(
            max_tries=self.execution.failure_policy.n_tries
        )
        self._subtask_container.load_data(all_enumerations)

        # execute Gold docking
        self._execute_gold()
