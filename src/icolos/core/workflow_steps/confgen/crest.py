import os
from typing import List

from pydantic import BaseModel
from rdkit import Chem
from copy import deepcopy

from icolos.utils.execute_external.crest import CrestExecutor

from icolos.utils.general.molecules import get_charge_for_molecule

from icolos.core.containers.compound import Enumeration, Conformer

from icolos.utils.enums.program_parameters import CrestEnum, CrestOutputEnum
from icolos.core.workflow_steps.step import _LE, _CTE
from icolos.core.workflow_steps.confgen.base import StepConfgenBase

_EE = CrestEnum()
_COE = CrestOutputEnum()


class StepCREST(StepConfgenBase, BaseModel):
    def __init__(self, **data):
        super().__init__(**data)

        # initialize the executor and test availability
        self._initialize_backend(executor=CrestExecutor)
        self._check_backend_availability()

    def _get_energies_from_XYZ(self, path) -> list:
        energies = []
        with open(path, "r") as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith(_COE.PREFIX_ENERGIES_XYZ):
                    energies.append(line.lstrip().rstrip())
        return energies

    def _parse_CREST_result(
        self, dir_path: str, enumeration: Enumeration
    ) -> List[Conformer]:
        """Function to parse the result from CREST."""
        # CREST will output a variety of files to "dir_path"
        conformers_sdf = os.path.join(dir_path, _COE.CREST_CONFORMERS_SDF)
        conformers_xyz = os.path.join(dir_path, _COE.CREST_CONFORMERS_XYZ)

        # as the energies are lost in the SDF output, we will add them as a tag
        energies = self._get_energies_from_XYZ(conformers_xyz)
        charge = str(
            get_charge_for_molecule(enumeration.get_molecule(), add_as_tag=False)
        )
        mol_supplier = Chem.SDMolSupplier(conformers_sdf, removeHs=False)
        result = []
        for mol_id, mol in enumerate(mol_supplier):
            mol.SetProp(_CTE.CONFORMER_ENERGY_TAG, energies[mol_id])
            mol.SetProp(_CTE.FORMAL_CHARGE_TAG, charge)
            result.append(Conformer(conformer=mol))
        return result

    def _set_formal_charge(self, parameters: dict, molecule: Chem.Mol) -> dict:
        charge = get_charge_for_molecule(molecule, add_as_tag=False)
        parameters[_EE.CREST_CHRG] = charge
        self._logger.log(f"Set charge for molecule to {charge}.", _LE.DEBUG)
        return parameters

    def _set_number_cores(self, parameters: dict) -> dict:
        """Function for parallelization of task, setting the number of cores to be used."""
        parameters[_EE.CREST_T] = int(self.execution.parallelization.jobs)
        return parameters

    def _prepare_settings(self, tmp_dir: str, enumeration: Enumeration) -> list:
        # first position is the input (SDF) file; the internal input at this stage is a molecule
        # -> write it to a temporary SDF file (undocumented input functionality) and add the path
        settings = [self._prepare_temp_input(tmp_dir, enumeration.get_molecule())]

        # add flags
        for flag in self.settings.arguments.flags:
            settings.append(flag)

        # add parameters
        parameters = deepcopy(self.settings.arguments.parameters)

        # update / over-write fields that need a specific value or are defined elsewhere
        parameters = self._set_number_cores(parameters)
        parameters = self._set_formal_charge(parameters, enumeration.get_molecule())

        # flatten the dictionary into a list for command-line execution
        for key in parameters.keys():
            settings.append(key)
            settings.append(parameters[key])
        return settings

    def execute(self):
        for compound in self.get_compounds():
            for enumeration in compound.get_enumerations():
                if not self._input_object_valid(enumeration):
                    continue

                # set up
                tmp_dir = self._move_to_temp_dir()

                # the call to CREST starts with the path to the input file, followed by arguments and flags
                settings = self._prepare_settings(tmp_dir, enumeration=enumeration)

                self._logger.log(
                    f"Executing CREST backend in folder {tmp_dir}.", _LE.DEBUG
                )
                result = self._backend_executor.execute(
                    command=_EE.CREST, arguments=settings, check=False
                )
                self._restore_working_dir()

                conformers = self._parse_CREST_result(tmp_dir, enumeration=enumeration)
                enumeration.clear_conformers()
                enumeration.add_conformers(conformers=conformers, auto_update=True)
                self._logger.log(
                    f"Executed CREST and obtained {len(conformers)} for enumeration {enumeration.get_index_string()}",
                    _LE.INFO,
                )

                self._remove_temporary(tmp_dir)
