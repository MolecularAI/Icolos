import os
import shutil
import tempfile
from typing import List, Tuple

from pydantic import BaseModel, Field
from rdkit import Chem
from copy import deepcopy

from icolos.utils.enums.step_enums import StepBaseEnum, StepGoldEnum
from icolos.utils.execute_external.autodockvina import AutoDockVinaExecutor
from icolos.utils.execute_external.gold import GoldExecutor
from icolos.utils.execute_external.openbabel import OpenBabelExecutor
from icolos.utils.general.icolos_exceptions import StepFailed
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
    protein_data: ConfigProteinData = Field(
        alias="PROTEIN DATA", default=None
    )


class GoldAdditional(BaseModel):
    configuration: GoldConfiguration = GoldConfiguration()
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
        try:
            result_tag_lines = conformer.GetProp(_ADE.REMARK_TAG).split("\n")
            result_line = [
                line for line in result_tag_lines if _ADE.RESULT_LINE_IDENTIFIER in line
            ][0]
            parts = result_line.split()
            docking_score = parts[_ADE.RESULT_LINE_POS_SCORE]
        except KeyError:
            return False
        conformer.SetProp(_SBE.ANNOTATION_TAG_DOCKING_SCORE, str(docking_score))
        return True

    def _generate_temporary_input_output_files(
        self, batch: List[List[Subtask]]
    ) -> Tuple[List[str], List[str], List[str], List[str]]:
        tmp_output_dirs = []
        tmp_input_paths = []
        tmp_output_paths = []
        enumeration_ids = []

        for next_subtask_list in batch:
            # for "AutoDock Vina", only single molecules can be handled so every sublist is
            # guaranteed at this stage to have only one element
            if len(next_subtask_list) > 1:
                self._logger.log(
                    f"Subtask list length for ADV is > 1 ({len(next_subtask_list)}), only the first element will be processed.",
                    _LE.WARNING,
                )
            subtask = next_subtask_list[0]

            # generate temporary input files and output directory
            cur_tmp_output_dir = tempfile.mkdtemp()
            _, cur_tmp_input_pdbqt = gen_tmp_file(
                suffix=".pdbqt", dir=cur_tmp_output_dir
            )
            _, cur_tmp_output_sdf = gen_tmp_file(suffix=".sdf", dir=cur_tmp_output_dir)

            # try to write the enumeration molecules out as PDBQT files
            enumeration = subtask.data
            mol = deepcopy(enumeration.get_molecule())
            if mol is None:
                shutil.rmtree(cur_tmp_output_dir)
                self._logger.log(
                    f"Enumeration {enumeration.get_index_string()} did not hold a valid RDkit molecule - skipped.",
                    _LE.DEBUG,
                )
                continue
            if not self._write_molecule_to_pdbqt(cur_tmp_input_pdbqt, mol):
                self._logger.log(
                    f"Could not generate PDBQT intermediate file from enumeration {enumeration.get_index_string()} - skipped.",
                    _LE.DEBUG,
                )
                continue

            # also store all the paths in case it succeeded -> these will be used later, failures will be ignored
            tmp_output_dirs.append(cur_tmp_output_dir)
            tmp_input_paths.append(cur_tmp_input_pdbqt)
            tmp_output_paths.append(cur_tmp_output_sdf)
            enumeration_ids.append(enumeration.get_index_string())

        return tmp_output_dirs, tmp_input_paths, tmp_output_paths, enumeration_ids

    def _execute_autodockvina(self):
        # get number of sublists in batch and initialize Parallelizer
        adv_parallelizer = Parallelizer(func=self._run_subjob)

        # continue until everything is successfully done or number of retries have been exceeded
        while self._subtask_container.done() is False:
            next_batch = self._get_sublists(get_first_n_lists=self._get_number_cores())

            # generate paths and initialize molecules (so that if they fail, this can be covered)
            (
                tmp_output_dirs,
                tmp_input_paths,
                tmp_output_paths,
                enumeration_ids,
            ) = self._generate_temporary_input_output_files(next_batch)

            # execute the current batch in parallel; hand over lists of parameters (will be handled by Parallelizer)
            # also increment the tries and set the status to "failed" (don't do that inside subprocess, as data is
            # copied, not shared!)
            _ = [sub.increment_tries() for element in next_batch for sub in element]
            _ = [sub.set_status_failed() for element in next_batch for sub in element]
            adv_parallelizer.execute_parallel(
                input_path_pdbqt=tmp_input_paths, output_path_sdf=tmp_output_paths
            )

            # parse the output of that particular batch and remove temporary files
            self._parse_adv_output_batch(
                tmp_input_paths=tmp_input_paths,
                tmp_output_paths=tmp_output_paths,
                enumeration_ids=enumeration_ids,
                next_batch=next_batch,
            )

            # clean-up
            self._remove_temporary(tmp_output_dirs)

            # print the progress for this execution
            self._log_execution_progress()

    def _parse_adv_output_batch(
        self,
        tmp_input_paths: List[str],
        tmp_output_paths: List[str],
        enumeration_ids: List[str],
        next_batch: List[List[Subtask]],
    ):

        for i in range(len(next_batch)):
            subtask = next_batch[i][0]
            tmp_output_path = tmp_output_paths[i]
            tmp_input_path = tmp_input_paths[i]
            enumeration_id = enumeration_ids[i]
            grid_id = self.adv_additional.grid_ids[0]
            grid_path = self.adv_additional.configuration.receptor_path

            # this is a protection against the case where empty (file size == 0 bytes) files are generated due to
            # a failure during docking
            if (
                not os.path.isfile(tmp_output_path)
                or os.path.getsize(tmp_output_path) == 0
            ):
                continue

            mol_supplier = Chem.SDMolSupplier(tmp_output_path, removeHs=False)
            for mol in mol_supplier:
                if mol is None:
                    continue
                cur_enumeration_name = str(mol.GetProp("_Name"))

                # add the information on the actual grid used
                mol.SetProp(_SBE.ANNOTATION_GRID_ID, str(grid_id))
                mol.SetProp(_SBE.ANNOTATION_GRID_PATH, str(grid_path))
                mol.SetProp(_SBE.ANNOTATION_GRID_FILENAME, os.path.basename(grid_path))

                # if no docking score is attached (i.e. the molecule is a receptor or so, skip it)
                if self._set_docking_score(mol) is not True:
                    continue

                # add molecule to the appropriate ligand
                for compound in self.get_compounds():
                    for enumeration in compound:
                        if enumeration.get_index_string() == enumeration_id:
                            new_conformer = Conformer(
                                conformer=mol,
                                conformer_id=None,
                                enumeration_object=enumeration,
                            )
                            enumeration.add_conformer(new_conformer, auto_update=True)
                            subtask.set_status_success()
                            break

    def _delay_file_system(self, path) -> bool:
        return self._wait_until_file_generation(
            path=path, interval_sec=2, maximum_sec=10
        )

    def _run_subjob(self, input_path_pdbqt: str, output_path_sdf: str):

        config = self.adv_additional.configuration

        # set up arguments list and execute
        _, tmp_pdbqt_docked = gen_tmp_file(
            suffix=".pdbqt", dir=os.path.dirname(input_path_pdbqt)
        )
        arguments = [
            _ADE.VINA_RECEPTOR,
            config.receptor_path,
            _ADE.VINA_LIGAND,
            input_path_pdbqt,
            _ADE.VINA_CPU,
            str(1),
            _ADE.VINA_SEED,
            config.seed,
            _ADE.VINA_OUT,
            tmp_pdbqt_docked,
            _ADE.VINA_CENTER_X,
            str(config.search_space.center_x),
            _ADE.VINA_CENTER_Y,
            str(config.search_space.center_y),
            _ADE.VINA_CENTER_Z,
            str(config.search_space.center_z),
            _ADE.VINA_SIZE_X,
            str(config.search_space.size_x),
            _ADE.VINA_SIZE_Y,
            str(config.search_space.size_y),
            _ADE.VINA_SIZE_Z,
            str(config.search_space.size_z),
            _ADE.VINA_NUM_MODES,
            config.number_poses,
        ]

        execution_result = self._backend_executor.execute(
            command=_ADE.VINA, arguments=arguments, check=True
        )
        self._delay_file_system(path=tmp_pdbqt_docked)

        # translate the parsed output PDBQT into an SDF
        arguments = [
            tmp_pdbqt_docked,
            _OBE.OBABEL_INPUTFORMAT_PDBQT,
            _OBE.OBABEL_OUTPUT_FORMAT_SDF,
            "".join([_OBE.OBABEL_O, output_path_sdf]),
        ]
        self._openbabel_executor.execute(
            command=_OBE.OBABEL, arguments=arguments, check=False
        )
        self._delay_file_system(path=output_path_sdf)

    def generate_config_file(self, path: str, ligand_files: List[str]):
        # very complicated custom format for GOLD; see file in "external_documentation/gold.conf" for details
        # TODO: constraints are not yet supported (and probably a few other fine-grained settings)
        def block_indent(block_name: str) -> str:
            return "".join([_SGE.BLOCK_INDENT, block_name])

        def empty_line() -> str:
            return ""

        def add_block(list_lines: List[str], block: dict):
            for key, value in block.items():
                if key == _SGE.LIGAND_DATA_FILE:
                    # this needs to be overwritten by the actual compound collection
                    self._logger.log(f"Do not set ligand data files explicitly, use the compound handover - skipping.",
                                     _LE.WARNING)
                    continue
                list_lines.append(' '.join([key, '=', str(value)]))

        config = [block_indent(_SGE.CONFIGURATION_START), empty_line()]

        # automatic settings
        config.append(block_indent(_SGE.AUTOMATIC_SETTINGS))
        add_block(config, self.gold_additional.configuration.automatic_settings.dict())
        config.append(empty_line())

        # population
        config.append(block_indent(_SGE.POPULATION))
        add_block(config, self.gold_additional.configuration.population.dict())
        config.append(empty_line())

        # genetic operators
        config.append(block_indent(_SGE.GENETIC_OPERATORS))
        add_block(config, self.gold_additional.configuration.genetic_operators.dict())
        config.append(empty_line())

        # flood fill
        config.append(block_indent(_SGE.FLOOD_FILL))
        add_block(config, self.gold_additional.configuration.flood_fill.dict())
        config.append(empty_line())

        # data files
        config.append(block_indent(_SGE.DATA_FILES))
        for ligand_file in ligand_files:
            config.append(' '.join([_SGE.LIGAND_DATA_FILE, ligand_file, "10"]))
        add_block(config, self.gold_additional.configuration.data_files.dict())
        config.append(empty_line())

        # flags
        config.append(block_indent(_SGE.FLAGS))
        add_block(config, self.gold_additional.configuration.flags.dict())
        config.append(empty_line())

        # termination
        config.append(block_indent(_SGE.TERMINATION))
        add_block(config, self.gold_additional.configuration.termination.dict())
        config.append(empty_line())

        # constraints
        config.append(block_indent(_SGE.CONSTRAINTS))
        add_block(config, self.gold_additional.configuration.constraints.dict())
        config.append(empty_line())

        # covalent bonding
        config.append(block_indent(_SGE.COVALENT_BONDING))
        add_block(config, self.gold_additional.configuration.covalent_bonding.dict())
        config.append(empty_line())

        # save options
        config.append(block_indent(_SGE.SAVE_OPTIONS))
        add_block(config, self.gold_additional.configuration.save_options.dict())
        config.append(empty_line())

        # fitness function settings
        config.append(block_indent(_SGE.FITNESS_FUNCTION_SETTINGS))
        add_block(config, self.gold_additional.configuration.fitness_function_settings.dict())
        config.append(empty_line())

        # protein data
        config.append(block_indent(_SGE.PROTEIN_DATA))
        add_block(config, self.gold_additional.configuration.protein_data.dict())
        config.append(empty_line())

        # write out
        with open(path, 'w') as f:
            for line in config:
                f.write(line + "\n")

    def execute(self):
        # Note: This step only supports one grid at a time, ensemble docking is taken care of at the workflow level!

        # in order to be able to efficiently execute ADV on the enumeration level, all of them have to be unrolled
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

        # execute ADV
        self._execute_autodockvina()


"""import os
import tempfile
import shutil
import multiprocessing
import pickle
from copy import deepcopy
from enum import Enum
from typing import Optional, List, Tuple, Dict, Any
from typing_extensions import Literal

import rdkit.Chem as Chem

from ccdc.docking import Docker as DockerGold
from ccdc.io import MoleculeReader, EntryWriter
from ccdc.protein import Protein
from pydantic import BaseModel

from dockstream.core.Schrodinger.Glide_docker import Parallelization
from dockstream.utils.enums.logging_enums import LoggingConfigEnum
from dockstream.utils.execute_external.Gold import GoldExecutor

from dockstream.core.docker import Docker
from dockstream.core.Gold.Gold_result_parser import GoldResultParser
from dockstream.utils.enums.Gold_enums import GoldLigandPreparationEnum
from dockstream.utils.enums.Gold_enums import GoldTargetKeywordEnum, GoldExecutablesEnum, GoldOutputEnum
from dockstream.utils.general_utils import gen_temp_file

from dockstream.utils.translations.molecule_translator import MoleculeTranslator
from dockstream.utils.dockstream_exceptions import DockingRunFailed


class GoldFitnessFunction(str, Enum):
    GOLDSCORE = "goldscore"
    CHEMSCORE = "chemscore"
    ASP = "asp"
    PLP = "plp"


class GoldResponseValue(str, Enum):
    FITNESS = "fitness"
    VALUE = "value"


class GoldParameters(BaseModel):
    prefix_execution: Optional[str] = None
    binary_location: Optional[str] = None
    receptor_paths: Optional[List[str]] = None
    time_limit_per_compound: Optional[int] = None
    parallelization: Optional[Parallelization]
    fitness_function: GoldFitnessFunction
    response_value: GoldResponseValue = GoldResponseValue.FITNESS
    early_termination: bool
    autoscale: float  # Autoscale percentage. very fast: 10, medium: 50, very slow: 100.
    ndocks: int = 10
    diverse_solutions: Optional[Tuple[bool, Optional[int], Optional[
        float]]] = None  # If diverse solutions is enabled this will be (True, cluster size, rmsd), otherwise (False, None, None). TODO: rework for GUI.

    def get(self, key: str) -> Any:
        return self.dict()[key]


_LP = GoldLigandPreparationEnum()
_TK = GoldTargetKeywordEnum()
_EE = GoldExecutablesEnum()
_ROE = GoldOutputEnum()
_LE = LoggingConfigEnum()


class Gold(Docker):

    backend: Literal["Gold"] = "Gold"
    parameters: GoldParameters

    _target_dict: Dict = None
    _Gold_executor: GoldExecutor = None
    _scoring_function_parameters: Dict[str, str] = None

    class Config:
        underscore_attrs_are_private = True

    def __init__(self, **run_parameters):
        # invoke base class's constructor first
        super().__init__(**run_parameters)

        # prepare and check Gold backend availability
        self._check_Gold_backend_availability()

        # parse the fitness function and response value set
        self._parse_fitness_function()

        # set the tag name for the scoring function and whether minimial or maximum values are better
        self._scoring_function_parameters = self._get_scoring_function_parameters()

    def _check_Gold_backend_availability(self):
        self._Gold_executor = GoldExecutor(
            prefix_execution=self.parameters.prefix_execution,
            binary_location=self.parameters.binary_location)
        if not self._Gold_executor.is_available():
            raise DockingRunFailed("Cannot initialize Gold docker, as Gold backend is not available - abort.")
        self._logger.log(f"Checked Gold backend availability (prefix_execution={self.parameters.prefix_execution}).",
                         _LE.DEBUG)

    def _parse_fitness_function(self):
        self._logger.log(
            f"Set fitness function to {self.parameters.fitness_function} and response value to {self.parameters.response_value}.",
            _LE.DEBUG)

    def _initialize_cavity(self, settings):
        # load the target dictionary specification and initialize the cavity
        target_path = self.parameters.receptor_paths[0]
        with open(target_path, "rb") as file:
            self._target_dict = pickle.load(file)
            self._logger.log(f"Loaded pickled cavity dictionary stored in file {target_path}.", _LE.DEBUG)
            if self._target_dict[_TK.VERSION] != _TK.CURRENT_VERSION:
                self._logger.log(
                    f"Version of pickled target ({self._target_dict[_TK.VERSION]}) is not the same as DockStream's ({_TK.CURRENT_VERSION}).",
                    _LE.WARNING)
        self._logger.log(f"Unpacked the target dictionary.", _LE.DEBUG)

        tmpdir = tempfile.mkdtemp()
        if self._target_dict[_TK.CAVITY_METHOD] == _TK.CAVITY_METHOD_REFERENCE:
            # write ligand to temporary file (ending copied over in settings)
            tmp_ref_ligand_path = gen_temp_file(suffix=self._target_dict[_TK.REFERENCE_LIGAND_FILENAME], dir=tmpdir)
            with open(tmp_ref_ligand_path, 'w') as file:
                for line in self._target_dict[_TK.REFERENCE_LIGAND]:
                    file.write(line)
                self._logger.log(
                    f"Wrote temporary ligand file {tmp_ref_ligand_path} with {len(self._target_dict[_TK.REFERENCE_LIGAND])} lines.",
                    _LE.DEBUG)

            # write target PDB to temporary file
            tmp_target_path = gen_temp_file(suffix=".pdb", dir=tmpdir)
            with open(tmp_target_path, 'w') as file:
                for line in self._target_dict[_TK.TARGET_PDB]:
                    file.write(line)
                self._logger.log(
                    f"Wrote temporary target file {tmp_target_path} with {len(self._target_dict[_TK.TARGET_PDB])} lines.",
                    _LE.DEBUG)

            # build the cavity
            ref_ligand = MoleculeReader(filename=tmp_ref_ligand_path)[0]
            self._prepare_protein(settings, tmp_target_path)
            protein = settings.proteins[0]
            settings.binding_site = settings.BindingSiteFromLigand(protein,
                                                                   ref_ligand,
                                                                   distance=self._target_dict[
                                                                       _TK.CAVITY_REFERENCE_DISTANCE])
            settings.reference_ligand_file = tmp_ref_ligand_path
        elif self._target_dict[_TK.CAVITY_METHOD] == _TK.CAVITY_METHOD_POINT:
            raise NotImplementedError
            # origin (x,x,x)
            # distance x
        else:
            raise DockingRunFailed("Specified cavity determination method not defined for GOLD.")
        self._logger.log(f"Initialized GOLD Protein.BindingSite with method {self._target_dict[_TK.CAVITY_METHOD]}.",
                         _LE.DEBUG)

    def add_molecules(self, molecules: list):

        mol_trans = MoleculeTranslator(self.ligands, force_mol_type=_LP.TYPE_RDKIT)
        mol_trans.add_molecules(molecules)
        self.ligands = mol_trans.get_as_rdkit()
        self._docking_performed = False

    def _generate_temporary_input_output_files(self, start_indices, sublists):
        # in case singletons are handed over, wrap them in a list for "zipping" later
        if not isinstance(start_indices, list):
            start_indices = [start_indices]
        if not isinstance(sublists, list):
            sublists = [sublists]

        tmp_output_dirs = []
        tmp_input_sdf_paths = []
        tmp_output_sdf_paths = []
        for start_index, sublist in zip(start_indices, sublists):
            # generate temporary input files and output directory
            cur_tmp_output_dir = tempfile.mkdtemp()
            cur_tmp_sdf = gen_temp_file(prefix=str(start_index), suffix=".sdf", dir=cur_tmp_output_dir)

            # write-out the temporary input file
            writer = Chem.SDWriter(cur_tmp_sdf)
            one_written = False
            for ligand in sublist:
                # initialize all ligands (as they could have failed)
                if ligand.get_molecule() is not None:
                    mol = deepcopy(ligand.get_molecule())
                    mol.SetProp("_Name", ligand.get_identifier())
                    one_written = True
                    writer.write(mol)
            writer.close()
            if one_written is False:
                if os.path.isdir(cur_tmp_output_dir):
                    shutil.rmtree(cur_tmp_output_dir)
                continue

            # add the path to which "_dock_subjob()" will write the result SDF
            output_sdf_path = gen_temp_file(prefix=str(start_index), suffix="_result.sdf", dir=cur_tmp_output_dir)
            tmp_output_dirs.append(cur_tmp_output_dir)
            tmp_output_sdf_paths.append(output_sdf_path)
            tmp_input_sdf_paths.append(cur_tmp_sdf)
        return tmp_output_dirs, tmp_input_sdf_paths, tmp_output_sdf_paths

    def _dock(self, number_cores: int):
        # partition ligands into sublists and distribute to processor cores for docking
        start_indices, sublists = self.get_sublists_for_docking(number_cores=number_cores)
        number_sublists = len(sublists)
        self._logger.log(f"Split ligands into {number_sublists} sublists for docking.", _LE.DEBUG)
        sublists_submitted = 0
        slices_per_iteration = min(number_cores, number_sublists)

        while sublists_submitted < len(sublists):
            upper_bound_slice = min((sublists_submitted + slices_per_iteration), len(sublists))
            cur_slice_start_indices = start_indices[sublists_submitted:upper_bound_slice]
            cur_slice_sublists = sublists[sublists_submitted:upper_bound_slice]

            # generate paths and initialize molecules (so that if they fail, this can be covered)
            tmp_output_dirs, tmp_input_sdf_paths, \
            tmp_output_sdf_paths = self._generate_temporary_input_output_files(cur_slice_start_indices,
                                                                               cur_slice_sublists)

            # run in parallel; wait for all subjobs to finish before proceeding
            processes = []
            for chunk_index in range(len(tmp_output_dirs)):
                p = multiprocessing.Process(target=self._dock_subjob, args=(tmp_input_sdf_paths[chunk_index],
                                                                            tmp_output_sdf_paths[chunk_index],
                                                                            tmp_output_dirs[chunk_index]))
                processes.append(p)
                p.start()
            for p in processes:
                p.join()

            # add the number of input sublists rather than the output temporary folders to account for cases where
            # entire sublists failed to produce an input structure
            sublists_submitted += len(cur_slice_sublists)

            # load the chunks and recombine the result; add conformations
            for chunk_index in range(len(tmp_output_dirs)):
                # this is a protection against the case where empty (file size == 0 bytes) files are generated due to
                # a failure during docking
                if not os.path.isfile(tmp_output_sdf_paths[chunk_index]) or os.path.getsize(
                        tmp_output_sdf_paths[chunk_index]) == 0:
                    continue

                for molecule in Chem.SDMolSupplier(tmp_output_sdf_paths[chunk_index], removeHs=False):

                    # it can happen, that ligands have "impossible chemistry" and will be loaded by RDkit as "None"
                    if molecule is None:
                        continue

                    # parse the molecule name (sorted by FITNESS not the score) which looks like:
                    # "0:0|0xa6enezm|sdf|1|dock6"
                    cur_conformer_name = str(molecule.GetProp("_Name")).split(sep='|')[0]

                    # add molecule to the appropriate ligand
                    for ligand in self.ligands:
                        if ligand.get_identifier() == cur_conformer_name:
                            ligand.add_conformer(molecule)
                            break

            # clean-up
            for path in tmp_output_dirs:
                shutil.rmtree(path)
            self._log_docking_progress(number_done=sublists_submitted, number_total=number_sublists)

        # update conformer names to contain the conformer id
        # -> <ligand_number>:<enumeration>:<conformer_number>
        reverse = True if self._get_scoring_function_parameters()[_ROE.BEST] == "max" else False
        for ligand in self.ligands:
            ligand.set_conformers(sorted(ligand.get_conformers(),
                                         key=lambda x: self._get_score_from_conformer(conformer=x),
                                         reverse=reverse))
            ligand.add_tags_to_conformers()

        # log any docking fails
        self._docking_fail_check()

        # generate docking results as dataframe
        result_parser = GoldResultParser(ligands=[ligand.get_clone() for ligand in self.ligands],
                                         fitness_function=self.parameters.fitness_function,
                                         response_value=self.parameters.response_value)
        self._df_results = result_parser.as_dataframe()

        # set docking flag
        self._docking_performed = True

    def _dock_subjob(self, sdf_ligand_path, path_sdf_results, tmp_output_dir):
        # 1) prepare Gold docker: (i) "clone" the docker instance, (ii) set remaining, ligang-specific settings and
        #                         (iii) initialize this chunk's ligands
        cur_docker = DockerGold()
        settings = cur_docker.settings
        settings.output_directory = tmp_output_dir
        settings.output_file = os.path.basename(path_sdf_results)
        settings.output_format = "sdf"
        settings.fitness_function = self.parameters.fitness_function
        settings.early_termination = self.parameters.early_termination
        settings.autoscale = self.parameters.autoscale

        if self.parameters.diverse_solutions is not None:
            settings.diverse_solutions = self.parameters.diverse_solutions

        self._initialize_cavity(settings)

        settings.add_ligand_file(sdf_ligand_path, ndocks=self.parameters.ndocks)

        # 2) write settings file
        settings_file_path = os.path.join(tmp_output_dir, _EE.GOLD_AUTO_CONFIG_NAME)
        settings.write(settings_file_path)
        with open(settings_file_path, 'r') as file:
            self._logger.log(f"Contents of configurations file {settings_file_path}:", _LE.DEBUG)
            for line in file:
                self._logger_blank.log(line.rstrip("\n"), _LE.DEBUG)

        # 3) run Gold docker
        arguments = [settings_file_path]
        execution_result = self._Gold_executor.execute(command=_EE.GOLD_AUTO,
                                                       arguments=arguments,
                                                       check=False)
        self._delay4file_system(path=path_sdf_results)
        self._logger.log(
            f"Finished sublist (input: {sdf_ligand_path}, output directory: {tmp_output_dir}), with return code '{execution_result.returncode}'.",
            _LE.DEBUG)

    def _prepare_protein(self, settings, tmp_protein_path):
        protein = Protein.from_file(tmp_protein_path)
        protein.remove_all_waters()
        protein.remove_unknown_atoms()
        protein.add_hydrogens()

        ligands = protein.ligands
        for l in ligands:
            protein.remove_ligand(l.identifier)
        protein_file_name = os.path.join(settings.output_directory, 'clean_%s.mol2' % protein.identifier)

        with EntryWriter(protein_file_name) as writer:
            writer.write(protein)
        settings.add_protein_file(protein_file_name)
        return ligands

    def write_docked_ligands(self, path, mode="all"):
        self._write_docked_ligands(path, mode, mol_type=_LP.TYPE_RDKIT)

    def _get_scoring_function_parameters(self):
        # get the appropriate name of the tag and whether minimal or maximal values are best for
        # the specified scoring function
        if self.parameters.response_value == GoldResponseValue.FITNESS:
            scoring_function_parameters = _ROE.DICT_FITNESS[self.parameters.fitness_function]
        elif self.parameters.response_value == GoldResponseValue.VALUE:
            scoring_function_parameters = _ROE.DICT_VALUE[self.parameters.fitness_function]
        else:
            raise ValueError("Parameter response value must be either fitness or value.")
        self._logger.log(f"Set scoring_function_parameters to {scoring_function_parameters} for obtaining the scores.",
                         _LE.DEBUG)
        return scoring_function_parameters

    def get_scores(self, best_only):

        return self._get_scores(best_only=best_only, best=self._scoring_function_parameters[_ROE.BEST])

    def write_result(self, path, mode="all"):

        return self._write_result(path=path, mode=mode, best=self._scoring_function_parameters[_ROE.BEST])

    def _get_score_from_conformer(self, conformer):
        return float(conformer.GetProp(self._scoring_function_parameters[_ROE.TAG]))

    def _sort_conformers(self, conformers: list, best=None) -> list:
        return super()._sort_conformers(conformers=conformers,
                                        best=self._scoring_function_parameters[_ROE.BEST])
"""
