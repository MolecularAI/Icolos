import os
from typing import List

from pydantic import BaseModel
from rdkit import Chem
from copy import deepcopy
from icolos.utils.execute_external.omega import OMEGAExecutor
from icolos.core.workflow_steps.step import _LE, _CTE
from icolos.utils.general.molecules import get_charge_for_molecule

from icolos.core.containers.compound import Enumeration, Conformer

from icolos.utils.enums.program_parameters import OMEGAEnum, OMEGAOutputEnum
from icolos.core.workflow_steps.confgen.base import StepConfgenBase

_EE = OMEGAEnum()
_COE = OMEGAOutputEnum()


class StepOmega(StepConfgenBase, BaseModel):
    def __init__(self, **data):
        super().__init__(**data)

        # initialize the executor and test availability
        self._initialize_backend(executor=OMEGAExecutor)
        self._check_backend_availability()

    def _parse_OMEGA_result(
        self, dir_path: str, enumeration: Enumeration
    ) -> List[Conformer]:
        # OMEGA will output a variety of files to "dir_path"
        conformers_sdf = os.path.join(dir_path, _COE.OUTPUT_SDF_NAME)

        # energies are added as a tag in the output
        mol_supplier = Chem.SDMolSupplier(conformers_sdf, removeHs=False)
        charge = str(
            get_charge_for_molecule(enumeration.get_molecule(), add_as_tag=False)
        )
        result = []
        for mol_id, mol in enumerate(mol_supplier):
            mol.SetProp(
                _CTE.CONFORMER_ENERGY_TAG, mol.GetProp(_COE.CLASSIC_ENERGY_OUTPUT_TAG)
            )
            mol.ClearProp(_COE.CLASSIC_ENERGY_OUTPUT_TAG)
            mol.SetProp(_CTE.FORMAL_CHARGE_TAG, charge)
            conf = Conformer(conformer=mol)
            result.append(conf)
        return result

    def _set_input_output_paths(self, parameters: dict, input_path: str) -> dict:
        # this is handled this way to overwrite any specifications from the user for the input / output paths as well
        parameters[_EE.CLASSIC_INPUT] = input_path
        parameters[_EE.CLASSIC_OUTPUT] = _COE.OUTPUT_SDF_NAME
        return parameters

    def _prepare_settings(self, tmp_dir: str, enumeration: Enumeration) -> list:
        # the first argument is the mode of binary "oeomega" (for now defaults to "classic")
        settings = [_EE.OMEGA_MODE_CLASSIC]

        # add flags
        # make sure, the energy tag is set as well
        for flag in self.settings.arguments.flags:
            settings.append(flag)
        if _EE.CLASSIC_SDENERGY not in settings:
            settings.append(_EE.CLASSIC_SDENERGY)

        # add parameters
        parameters = deepcopy(self.settings.arguments.parameters)

        # update / over-write fields that need a specific value or are defined elsewhere
        parameters = self._set_input_output_paths(
            parameters=parameters,
            input_path=self._prepare_temp_input(tmp_dir, enumeration.get_molecule()),
        )

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
                settings = self._prepare_settings(tmp_dir, enumeration=enumeration)

                # execution
                self._logger.log(
                    f"Executing OMEGA backend in folder {tmp_dir}.", _LE.DEBUG
                )
                result = self._backend_executor.execute(
                    command=_EE.OMEGA, arguments=settings, check=False
                )
                self._restore_working_dir()

                # parsing
                conformers = self._parse_OMEGA_result(tmp_dir, enumeration=enumeration)
                enumeration.clear_conformers()
                enumeration.add_conformers(conformers=conformers, auto_update=True)
                self._logger.log(
                    f"Completed OMEGA for enumeration {enumeration.get_index_string()}, added {len(conformers)} conformers.",
                    _LE.INFO,
                )

                # clean-up
                self._remove_temporary(tmp_dir)
