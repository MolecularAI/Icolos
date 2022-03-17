import os
from tempfile import mkdtemp

from pydantic import BaseModel
from rdkit import Chem
from copy import deepcopy
from typing import List, Tuple
from icolos.utils.execute_external.xtb import XTBExecutor

from icolos.utils.general.molecules import get_charge_for_molecule

from icolos.core.containers.compound import Conformer

from icolos.utils.enums.program_parameters import XTBEnum, XTBOutputEnum
from icolos.core.workflow_steps.step import _LE, _CTE
from icolos.core.workflow_steps.confgen.base import StepConfgenBase
from icolos.utils.general.parallelization import SubtaskContainer, Parallelizer

_EE = XTBEnum()
_COE = XTBOutputEnum()


class StepXTB(StepConfgenBase, BaseModel):
    def __init__(self, **data):
        super().__init__(**data)

        # initialize the executor and test availability
        self._initialize_backend(executor=XTBExecutor)
        self._check_backend_availability()

    def _parse_XTB_result(self, tmp_dirs: List, conformers: List[Conformer]):
        # XTB will output a variety of files to "dir_path"
        results = []
        for dir_path, conformer in zip(tmp_dirs, conformers):
            optimized_conformer_sdf = os.path.join(dir_path, _COE.XTBOPT_SDF)
            enum = conformer.get_enumeration_object()

            # as the energies are added as a tag, but we will use ours
            # note, that XTB is called to operate on one conformer at a time (which we will return here)
            mol_supplier = Chem.SDMolSupplier(optimized_conformer_sdf, removeHs=False)
            try:
                for mol in mol_supplier:
                    mol.SetProp(
                        _CTE.CONFORMER_ENERGY_TAG, mol.GetProp(_COE.TOTAL_ENERGY_TAG)
                    )
                    mol.ClearProp(_COE.TOTAL_ENERGY_TAG)
                    mol.SetProp(
                        _CTE.FORMAL_CHARGE_TAG, str(get_charge_for_molecule(mol))
                    )
                    enum.add_conformer(Conformer(conformer=mol), auto_update=True)
                results.append(_COE.SUCCESS)
            except:
                self._logger.log(
                    f"Failed to parse XTB results for conformer {conformer.get_index_string()}.",
                    _LE.WARNING,
                )
                results.append(_COE.FAILURE)
        return results

    def _prepare_batch(self, batch) -> Tuple:
        # first position is the input (SDF) file; the internal input at this stage is a list of molecules
        # -> write it to a temporary SDF file (undocumented input functionality) and add the path

        tmp_dirs = []
        input_files = []
        charges = []
        conformers = []
        for next_subtask_list in batch:
            tmp_dir = mkdtemp()
            tmp_dirs.append(tmp_dir)
            for (
                subtask
            ) in (
                next_subtask_list
            ):  # enforced as one since xtb can't handle multiple files in one call
                conformer = subtask.data
                conformers.append(conformer)
                input_file = self._prepare_temp_input(tmp_dir, conformer.get_molecule())
                charge = get_charge_for_molecule(conformer.get_molecule())

                charges.append(charge)
                input_files.append(input_file)
        return tmp_dirs, input_files, charges, conformers

    def _prepare_arguments(self, settings: List) -> List:

        # add flags
        for flag in self.settings.arguments.flags:
            settings.append(flag)

        # add parameters
        parameters = deepcopy(self.settings.arguments.parameters)

        # flatten the dictionary into a list for command-line execution
        for key in parameters.keys():
            settings.append(key)
            settings.append(parameters[key])
        return settings

    def _run_subjob(self, tmp_dir: str, input_file: str, charge: int) -> None:

        work_dir = os.getcwd()
        os.chdir(tmp_dir)

        arguments = [input_file, _EE.XTB_P, charge]
        arguments = self._prepare_arguments(
            arguments
        )  # add additional parameters from config

        result = self._backend_executor.execute(
            command=_EE.XTB, arguments=arguments, check=False
        )
        # for line in result.stdout.split("\n"):
        #     self._logger_blank.log(line, _LE.DEBUG)
        #     # print(line)
        os.chdir(work_dir)

    def _execute_xtb(self):
        xtb_parallelizer = Parallelizer(func=self._run_subjob)
        n = 1

        tmp_dirs = None
        while self._subtask_container.done() is False:

            next_batch = self._get_sublists(get_first_n_lists=self._get_number_cores())
            tmp_dirs, input_files, charges, conformers = self._prepare_batch(next_batch)

            _ = [sub.increment_tries() for element in next_batch for sub in element]
            _ = [sub.set_status_failed() for element in next_batch for sub in element]

            self._logger.log(f"Executing xtb for batch {n}.", _LE.DEBUG)

            xtb_parallelizer.execute_parallel(
                tmp_dir=tmp_dirs,
                input_file=input_files,
                charge=charges,
            )

            results = self._parse_XTB_result(tmp_dirs, conformers)
            for sublist, result in zip(next_batch, results):
                assert len(sublist) == 1
                # TODO: this only works if max length sublist == 1, fine for now as that is all turbomole can handle
                for task in sublist:
                    if result == _COE.SUCCESS:
                        task.set_status_success()
                    else:
                        task.set_status_failed()

            n += 1
        self._remove_temporary(tmp_dirs)

    def execute(self):
        all_conformers = []
        for compound in self.get_compounds():
            for enumeration in compound.get_enumerations():
                if enumeration.get_conformers():
                    for conformer in enumeration.get_conformers():
                        all_conformers.append(conformer)
                enumeration.clear_conformers()
        self.execution.parallelization.max_length_sublists = 1
        self._subtask_container = SubtaskContainer(
            max_tries=self.execution.failure_policy.n_tries
        )
        self._subtask_container.load_data(all_conformers)
        self._execute_xtb()
        self._logger.log(
            f"Completed execution of XTB for {len(all_conformers)} conformers.",
            _LE.INFO,
        )
