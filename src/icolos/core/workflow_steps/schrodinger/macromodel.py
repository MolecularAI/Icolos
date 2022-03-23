import os
import subprocess
from typing import Tuple, List

from pydantic import BaseModel, PrivateAttr
from rdkit import Chem

from icolos.core.workflow_steps.schrodinger.base import StepSchrodingerBase
from icolos.utils.execute_external.macromodel import MacromodelExecutor

from icolos.utils.general.molecules import get_charge_for_molecule

from icolos.core.containers.compound import Enumeration, Conformer

from icolos.utils.enums.program_parameters import (
    MacromodelEnum,
)
from icolos.utils.enums.step_enums import StepMacromodelEnum
from icolos.core.workflow_steps.step import _LE, _CTE
from icolos.core.step_utils.sdconvert_util import SDConvertUtil

_EE = MacromodelEnum()
_MMSE = StepMacromodelEnum()


class StepMacromodel(StepSchrodingerBase, BaseModel):
    class Config:
        underscore_attrs_are_private = True

    _sdconvert_util = PrivateAttr()

    def __init__(self, **data):
        super().__init__(**data)

        # initialize the executor and test availability
        self._initialize_backend(executor=MacromodelExecutor)
        self._check_backend_availability()

        # prepare sdconvert utility
        self._sdconvert_util = SDConvertUtil(
            prefix_execution=self.execution.prefix_execution,
            binary_location=self.execution.binary_location,
        )

        # extend parameters with the COM file default, if not present
        if _MMSE.COM_FILE not in self.settings.arguments.parameters.keys():
            self.settings.arguments.parameters[_MMSE.COM_FILE] = _MMSE.COM_FILE_DEFAULT

    def _execute_macromodel(self, com_file: str) -> subprocess.CompletedProcess:
        self._logger.log(
            f"Executing MacroModel backend for com_file {com_file}.", _LE.DEBUG
        )
        arguments = []
        for key in self.settings.arguments.parameters.keys():
            # TODO: disentangle "special behaviour" for this key - move the com_file specification to a separate block
            #       in the configuration
            if key != _MMSE.COM_FILE:
                arguments.append(key)
                arguments.append(str(self.settings.arguments.parameters[key]))
        for flag in self.settings.arguments.flags:
            arguments.append(str(flag))
        arguments.append(com_file)
        self._apply_token_guard()
        result = self._backend_executor.execute(
            command=_EE.MACROMODEL, arguments=arguments, check=True
        )
        return result

    def _set_formal_charge(self, parameters: dict, molecule: Chem.Mol) -> dict:
        charge = get_charge_for_molecule(molecule)
        parameters[_EE.XTB_CHRG] = charge
        self._logger.log(f"Set charge for molecule to {charge}.", _LE.DEBUG)
        return parameters

    def _prepare_file_paths(self, tmp_dir: str) -> Tuple[str, str, str]:
        # generate the paths to the temporary files
        mae_input = os.path.join(tmp_dir, _MMSE.MAE_INPUT)
        mae_output = os.path.join(tmp_dir, _MMSE.MAE_OUTPUT)
        sdf_output = os.path.join(tmp_dir, _MMSE.SDF_OUTPUT)

        return mae_input, mae_output, sdf_output

    def _prepare_settings_file(self, tmp_dir: str) -> str:
        path_settings_file = os.path.join(tmp_dir, _MMSE.COM_FILE_PATH)

        # join the input and output paths (at the beginning of the COM file) and the
        # settings from either the default or the configuration together
        complete_com = "\n".join(
            [
                os.path.join(tmp_dir, _MMSE.MAE_INPUT),
                os.path.join(tmp_dir, _MMSE.MAE_OUTPUT),
                self.settings.arguments.parameters[_MMSE.COM_FILE],
            ]
        )
        with open(path_settings_file, "w") as f:
            f.writelines(complete_com)
        return path_settings_file

    def _prepare_run_files(
        self, tmp_dir: str, enumeration: Enumeration
    ) -> Tuple[str, str, str, str, str]:
        # generate the file paths (NOT populated yet)
        mae_input, mae_output, sdf_output = self._prepare_file_paths(tmp_dir)

        # write the input SDF file and translate it into Schrodingers native MAE format
        sdf_input = self._prepare_temp_input(tmp_dir, enumeration.get_molecule())
        self._sdconvert_util.sdf2mae(sdf_input, mae_input)

        # write out the settings file
        com_file = self._prepare_settings_file(tmp_dir)

        return sdf_input, mae_input, mae_output, sdf_output, com_file

    def _parse_macromodel_result(
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

    def execute(self):
        for compound in self.get_compounds():
            for enumeration in compound.get_enumerations():
                if not self._input_object_valid(enumeration):
                    continue

                # set up
                tmp_dir = self._move_to_temp_dir()

                # get the paths to the MAE and SDF input and output files and the COM file (settings)
                (
                    sdf_input,
                    mae_input,
                    mae_output,
                    sdf_output,
                    com_file,
                ) = self._prepare_run_files(tmp_dir=tmp_dir, enumeration=enumeration)

                # execute MacroModel, obtain the output SDF and switch back the working directory to what it was before
                _ = self._execute_macromodel(com_file=com_file)
                os.listdir(tmp_dir)
                self._sdconvert_util.mae2sdf(mae_file=mae_output, sdf_file=sdf_output)
                self._restore_working_dir()

                # parse output
                conformers = self._parse_macromodel_result(sdf_output, enumeration)
                enumeration.clear_conformers()
                enumeration.add_conformers(conformers=conformers, auto_update=True)
                self._logger.log(
                    f"Executed MacroModel and obtained {len(conformers)} conformers for enumeration {enumeration.get_index_string()}.",
                    _LE.INFO,
                )

                self._remove_temporary(tmp_dir)
