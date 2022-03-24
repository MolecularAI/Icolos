import os
import tempfile
from typing import List

from pydantic import BaseModel
from rdkit import Chem

from icolos.core.workflow_steps.schrodinger.base import StepSchrodingerBase
from icolos.utils.execute_external.ligprep import LigprepExecutor
from icolos.utils.general.files_paths import gen_tmp_file

from icolos.utils.general.molecules import get_charge_for_molecule
from icolos.core.containers.compound import Enumeration, Conformer, get_compound_by_id

from icolos.utils.enums.program_parameters import (
    LigprepEnum,
)
from icolos.utils.enums.step_enums import StepLigprepEnum
from icolos.core.workflow_steps.step import _LE, _CTE
from icolos.utils.general.parallelization import SubtaskContainer, Parallelizer, Subtask
from icolos.utils.general.print_log import print_log_file
from icolos.utils.smiles import to_smiles

_EE = LigprepEnum()
_SLE = StepLigprepEnum()


class StepLigprep(StepSchrodingerBase, BaseModel):
    """
    Interface to the LigPrep binary for ligand embedding
    """

    class Config:
        underscore_attrs_are_private = True

    def __init__(self, **data):
        super().__init__(**data)

        # initialize the executor and test availability
        self._initialize_backend(executor=LigprepExecutor)
        self._check_backend_availability()

    def _prepare_ligprep_arguments(self) -> list:
        arguments_list = []

        # add user-specified command-line settings (if provided); note, that empty dictionaries evaluate
        # to False
        if len(self.settings.arguments.flags) > 0:
            for flag in self.settings.arguments.flags:
                arguments_list.append(str(flag))
        if self.settings.arguments.parameters:
            for key in self.settings.arguments.parameters.keys():
                if key == _EE.LIGPREP_F:
                    self._logger.log(
                        'Removing "-f" parameter for Ligprep arguments - filter file settings need to be specified in the "additional" block directly.',
                        _LE.WARNING,
                    )
                    continue

                arguments_list.append(key)
                if (
                    self.settings.arguments.parameters[key] is not None
                    and self.settings.arguments.parameters[key] != ""
                ):
                    arguments_list.append(str(self.settings.arguments.parameters[key]))

        # add default settings, that are not exposed to the user yet
        if _EE.LIGPREP_HOST not in arguments_list:
            arguments_list.append(_EE.LIGPREP_HOST)
            arguments_list.append(_EE.LIGPREP_HOST_LOCALHOST)
        arguments_list.append(_EE.LIGPREP_WAIT)
        arguments_list = arguments_list + [_EE.LIGPREP_NJOBS, 1]

        return arguments_list

    def _generate_temporary_input_output_files(self, batch: List[List[Subtask]]):
        tmp_output_dirs = []
        tmp_input_smi_paths = []
        tmp_input_filter_paths = []
        tmp_output_sdf_paths = []
        dict_original_smiles = {}

        for next_subtask_list in batch:
            # generate temporary input files and output directory
            cur_tmp_output_dir = tempfile.mkdtemp()
            _, cur_tmp_smi = gen_tmp_file(suffix=".smi", dir=cur_tmp_output_dir)
            _, cur_tmp_filter = gen_tmp_file(suffix=".lff", dir=cur_tmp_output_dir)

            # write smiles to temporary file as "Ligprep" backend
            with open(cur_tmp_smi, "w") as f:
                for subtask in next_subtask_list:
                    enumeration = subtask.data
                    dict_original_smiles[
                        enumeration.get_index_string()
                    ] = enumeration.get_original_smile()
                    f.write(
                        enumeration.get_original_smile()
                        + " "
                        + enumeration.get_index_string()
                        + "\n"
                    )

            # add the path to which "_dock_subjob()" will write the result SDF
            _, output_sdf_path = gen_tmp_file(
                suffix="_result.sdf", dir=cur_tmp_output_dir
            )

            # add the temporary paths
            tmp_output_dirs.append(cur_tmp_output_dir)
            tmp_input_smi_paths.append(cur_tmp_smi)
            tmp_input_filter_paths.append(cur_tmp_filter)
            tmp_output_sdf_paths.append(output_sdf_path)
        return (
            tmp_output_dirs,
            tmp_input_smi_paths,
            tmp_output_sdf_paths,
            tmp_input_filter_paths,
            dict_original_smiles,
        )

    def _execute_ligprep(self):
        ligprep_parallelizer = Parallelizer(func=self._run_subjob)

        # continue until everything is successfully done or number of retries have been exceeded
        while self._subtask_container.done() is False:
            next_batch = self._get_sublists(get_first_n_lists=self._get_number_cores())

            # generate paths and initialize molecules (so that if they fail, this can be covered)
            (
                tmp_output_dirs,
                tmp_input_smi_paths,
                tmp_output_sdf_paths,
                tmp_input_filter_paths,
                dict_original_smiles,
            ) = self._generate_temporary_input_output_files(next_batch)

            # call "token guard" method (only executed, if block is specified in the configuration), which will wait
            # with the execution if not enough tokens are available at the moment
            self._apply_token_guard()

            # execute the current batch in parallel; hand over lists of parameters (will be handled by Parallelizer)
            # also increment the tries and set the status to "failed" (don't do that inside subprocess, as data is
            # copied, not shared!)
            _ = [sub.increment_tries() for element in next_batch for sub in element]
            _ = [sub.set_status_failed() for element in next_batch for sub in element]
            ligprep_parallelizer.execute_parallel(
                smi_ligand_path=tmp_input_smi_paths,
                path_sdf_results=tmp_output_sdf_paths,
                tmp_output_dir=tmp_output_dirs,
                tmp_input_filter=tmp_input_filter_paths,
                sublist=next_batch,
            )

            # parse the output of that particular batch and remove temporary files
            self._parse_ligprep_output(
                tmp_output_sdf_paths, dict_original_smiles, next_batch
            )
            self._remove_temporary(tmp_output_dirs)

            # print the progress for this execution
            self._log_execution_progress()

    def _parse_ligprep_output(
        self,
        tmp_output_sdf_paths: List[str],
        dict_original_smiles: dict,
        batch: List[List[Subtask]],
    ):
        # TODO: refactor that
        def _update_subtask(sublist: List[Subtask], enum_identifier: str):
            for task in sublist:
                if task.data.get_index_string() == enum_identifier:
                    task.set_status_success()

        for i in range(len(tmp_output_sdf_paths)):
            # get input and output paths and check the files are there
            path_sdf_results = tmp_output_sdf_paths[i]
            cur_sublist = batch[i]
            if (
                not os.path.isfile(path_sdf_results)
                or os.path.getsize(path_sdf_results) == 0
            ):
                continue

            mol_supplier = Chem.SDMolSupplier(path_sdf_results, removeHs=False)
            for mol in mol_supplier:
                # Ligprep adds a "-1" to "-[N]" to the names in the variants tag; this tag is always added
                # alternatively, the "_Name" property could be loaded
                # TODO: add loading only the most likely tautomer here (based on _SLE.LIGPREP_TAUTOMER_PROBABILITY)
                if mol is not None and mol.HasProp(_SLE.LIGPREP_VARIANTS):
                    identifier, _ = mol.GetProp(_SLE.LIGPREP_VARIANTS).split("-")
                    compound_id, enumeration_id = identifier.split(":")
                    compound = get_compound_by_id(
                        self.get_compounds(), int(compound_id)
                    )
                    enumeration = Enumeration(
                        compound_object=compound,
                        smile=to_smiles(mol),
                        original_smile=dict_original_smiles[identifier],
                        molecule=mol,
                    )
                    compound.add_enumeration(enumeration, auto_update=True)
                    _update_subtask(cur_sublist, enum_identifier=identifier)
                else:
                    self._logger.log(
                        f"Skipped molecule when loading as specified property {_SLE.LIGPREP_VARIANTS} could not be found - typically, this indicates that ligprep could not embed the molecule.",
                        _LE.WARNING,
                    )

    def _add_filtering(self, arguments: list, tmp_input_filter: str) -> list:
        filter_file_settings = self.settings.additional.get(_SLE.FILTER_FILE, None)
        if filter_file_settings is not None:
            filter_file = open(tmp_input_filter, "w")
            for key in filter_file_settings.keys():
                filter_file.write(
                    f"{key}                             {filter_file_settings[key]}\n"
                )
            filter_file.close()
            arguments = arguments + [_EE.LIGPREP_F, tmp_input_filter]
        return arguments

    def _run_subjob(
        self,
        smi_ligand_path: str,
        path_sdf_results: str,
        tmp_output_dir: str,
        tmp_input_filter: str,
        sublist: List[Subtask],
    ):
        # 1) increase the sublist "tries" and set status to "failed"
        _ = [task.increment_tries() for task in sublist]
        _ = [task.set_status_failed() for task in sublist]

        # 2) change to directory, to be able to use relative paths (to compensate for Schrodinger bug with AWS)
        working_dir = os.getcwd()
        os.chdir(tmp_output_dir)

        # 3) prepare "Ligprep" arguments
        arguments = self._prepare_ligprep_arguments()
        arguments = self._add_filtering(arguments, tmp_input_filter)
        arguments = arguments + [
            _EE.LIGPREP_INPUT_ISMI,
            os.path.basename(smi_ligand_path),
        ]
        arguments = arguments + [
            _EE.LIGPREP_OUTPUT_OSD,
            os.path.basename(path_sdf_results),
        ]

        # 4) run "Ligprep" backend and add log file to "debug" mode logging
        result = self._backend_executor.execute(
            command=_EE.LIGPREP,
            arguments=arguments,
            location=tmp_output_dir,
            check=False,
        )

        self._logger.log(
            f"Executed Ligprep backend (output file: {path_sdf_results}).", _LE.DEBUG
        )
        path_tmp_log = os.path.join(
            tmp_output_dir,
            "".join(
                [
                    os.path.splitext(os.path.basename(path_sdf_results))[0],
                    _EE.LIGPREP_LOG_ENDING,
                ]
            ),
        )
        print_log_file(path=path_tmp_log, logger=self._logger, level=_LE.DEBUG)

        # 5) revert back to working directory
        os.chdir(working_dir)

    def _parse_ligprep_result(
        self, sdf_output: str, enumeration: Enumeration
    ) -> List[Conformer]:
        charge = str(
            get_charge_for_molecule(enumeration.get_molecule(), add_as_tag=False)
        )
        mol_supplier = Chem.SDMolSupplier(sdf_output, removeHs=False)
        conformers = []
        for mol_id, mol in enumerate(mol_supplier):
            # note, that formal charge information would be kept if available before (i.e. it retains tags)
            mol.SetProp(_CTE.FORMAL_CHARGE_TAG, charge)
            conformers.append(Conformer(conformer=mol))
        return conformers

    def _log_execution(self, initial_enum_number: int):
        number_enumerations_after = 0
        for compound in self.get_compounds():
            number_enumerations_after += len(compound.get_enumerations())
        self._logger.log(
            f"Executed LigPrep for {initial_enum_number} input enumerations, resulting in {number_enumerations_after} output enumerations.",
            _LE.INFO,
        )
        for compound in self.get_compounds():
            for enumeration in compound.get_enumerations():
                self._logger.log(
                    f"Added enumeration {enumeration.get_index_string()} with smile {enumeration.get_smile()}.",
                    _LE.DEBUG,
                )

    def execute(self):
        # in order to be able to efficiently execute Ligprep on the enumeration level, all of them have to be unrolled
        # Note: As they retain their respective Compound object, the attribution later on is simple
        all_enumerations = []
        for compound in self.get_compounds():
            all_enumerations = all_enumerations + compound.get_enumerations()
            compound.clear_enumerations()
        # TODO: we will use the "original_smile" of the enumeration to start the embedding; make sure it exists

        # split into sublists, according to the settings
        self._subtask_container = SubtaskContainer(
            max_tries=self.execution.failure_policy.n_tries
        )
        self._subtask_container.load_data(all_enumerations)

        self._execute_ligprep()
        self._log_execution(initial_enum_number=len(all_enumerations))
