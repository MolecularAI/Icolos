import os
import shutil
import tempfile
from typing import List, Tuple

from pydantic import BaseModel, Field
from rdkit import Chem
from copy import deepcopy

from icolos.utils.enums.step_enums import StepAutoDockVinaEnum, StepBaseEnum
from icolos.utils.execute_external.autodockvina import AutoDockVinaExecutor
from icolos.utils.execute_external.openbabel import OpenBabelExecutor
from icolos.utils.general.icolos_exceptions import StepFailed
from icolos.utils.general.files_paths import gen_tmp_file
from icolos.core.containers.compound import Conformer
from icolos.utils.enums.program_parameters import AutoDockVinaEnum, OpenBabelEnum
from icolos.core.workflow_steps.step import _LE, StepBase
from icolos.utils.general.parallelization import Subtask, SubtaskContainer, Parallelizer

_SBE = StepBaseEnum
_ADE = AutoDockVinaEnum()
_OBE = OpenBabelEnum()
_SAE = StepAutoDockVinaEnum()


class ADVSearchSpace(BaseModel):
    center_x: float = Field(alias="--center_x", default=None)
    center_y: float = Field(alias="--center_y", default=None)
    center_z: float = Field(alias="--center_z", default=None)
    size_x: float = Field(alias="--size_x", default=15.0)
    size_y: float = Field(alias="--size_y", default=15.0)
    size_z: float = Field(alias="--size_z", default=15.0)


class ADVConfiguration(BaseModel):
    seed: int = 42
    number_poses: int = 1
    search_space: ADVSearchSpace = ADVSearchSpace()
    receptor_path: str = None


class ADVAdditional(BaseModel):
    configuration: ADVConfiguration = ADVConfiguration()
    grid_ids: List[str] = ["grid0"]


class StepAutoDockVina(StepBase, BaseModel):

    _openbabel_executor: OpenBabelExecutor = None
    adv_additional: ADVAdditional = None

    def __init__(self, **data):
        super().__init__(**data)

        # initialize the executor and test availability
        self._initialize_backend(executor=AutoDockVinaExecutor)
        self._check_backend_availability()

        # initialize the executor for all "OpenBabel"
        self._openbabel_executor = OpenBabelExecutor()
        if not self._openbabel_executor.is_available():
            raise StepFailed(
                "AutoDock Vina requires OpenBabel execution, initialization failed."
            )

        # set ADV specific settings and ensure that each molecule gets its own sublist
        self.adv_additional = ADVAdditional(**self.settings.additional)
        self.execution.parallelization.max_length_sublists = 1

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

    def _write_molecule_to_pdbqt(self, path: str, molecule: Chem.Mol) -> bool:
        # generate temporary copy as PDB
        _, tmp_pdb = gen_tmp_file(suffix=".pdb", dir=os.path.dirname(path))
        Chem.MolToPDBFile(molecule, filename=tmp_pdb)

        # translate the pdb into a pdbqt including partial charges
        # Note: In contrast to the target preparation,
        # we will use a tree-based flexibility treatment here -
        # thus, the option "-xr" is NOT used.
        arguments = [
            tmp_pdb,
            _OBE.OBABEL_OUTPUT_FORMAT_PDBQT,
            "".join([_OBE.OBABEL_O, path]),
            _OBE.OBABEL_PARTIALCHARGE,
            _OBE.OBABEL_PARTIALCHARGE_GASTEIGER,
        ]
        self._openbabel_executor.execute(
            command=_OBE.OBABEL, arguments=arguments, check=False
        )

        if os.path.exists(path):
            return True
        else:
            return False

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
