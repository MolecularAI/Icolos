from icolos.utils.enums.step_enums import StepGromacsEnum
from icolos.utils.enums.program_parameters import SchrodingerExecutablesEnum
from icolos.core.workflow_steps.schrodinger.base import StepSchrodingerBase
from icolos.utils.execute_external.execute import Executor
from pydantic import BaseModel
from copy import deepcopy
from icolos.core.workflow_steps.step import _LE
import os
import pandas as pd


class StepResidueScanning(StepSchrodingerBase, BaseModel):
    """
    Interface to Schrodinger's PrepWizard program for protein prep
    """

    def __init__(self, **data):
        super().__init__(**data)
        self._initialize_backend(executor=Executor)
        self._check_backend_availability()

    def _parse_args(self):
        parameters = deepcopy(self.settings.arguments.parameters)
        arguments = []
        if len(self.settings.arguments.flags) > 0:
            for flag in self.settings.arguments.flags:
                arguments.append(str(flag))
        if parameters:
            for key in parameters.keys():
                arguments.append(key)
                if parameters[key] is not None and parameters[key] != "":
                    arguments.append(str(parameters[key]))
        input_file = self.data.generic.get_argument_by_extension("maegz")
        arguments.append(input_file)
        print(arguments)
        return arguments

    def _parse_output(self, tmp_dir: str):
        # run the conversion script, generate csv for the mutants, attach to compounds
        out_file = [f for f in os.listdir(tmp_dir) if f.endswith("out.maegz")][0]
        command = f"$SCHRODINGER/utilities/proplister {out_file} -p s_bioluminate_Mutations -p r_bioluminate_delta_Stability -p r_bioluminate_delta_Affinity -c -o props_out.csv"
        self._backend_executor.execute(
            command, arguments=[], check=True, location=tmp_dir
        )

        # parse csv file, attach to compounds
        compounds = self.get_compounds()
        props_df = pd.read_csv("props_out.csv")
        for comp in compounds:
            mol = comp.get_enumerations()[0].get_molecule()

            # find the row by the name
            mol_name = mol.GetProp("_Name")
            row = props_df.loc[
                props_df["s_bioluminate_Mutations"].str.contains(mol_name)
            ]
            try:
                mol.SetProp(
                    "r_bioluminate_delta_Stability",
                    str(float(row["r_bioluminate_delta_Stability"])),
                )
                mol.SetProp(
                    "r_bioluminate_delta_Affinity",
                    str(float(row["r_bioluminate_delta_Affinity"])),
                )
                self._logger.log(
                    f"Parsed properties for mol {mol_name}: {mol.GetProp('r_bioluminate_delta_Affinity')}, {mol.GetProp('r_bioluminate_delta_Stability')}",
                    _LE.DEBUG,
                )
            except KeyError:
                self._logger.log("Failed to parse results df", _LE.WARNING)
                mol.SetProp(
                    "r_bioluminate_delta_Stability",
                    str(0.00),
                )
                mol.SetProp(
                    "r_bioluminate_delta_Affinity",
                    str(0.00),
                )

    def execute(self):
        tmp_dir = self._make_tmpdir()
        print(tmp_dir)
        args = self._parse_args()
        self.data.generic.write_out_all_files(tmp_dir)

        command = f"$SCHRODINGER/run residue_scanning_backend.py"
        self._logger.log("executing residue scanning", _LE.DEBUG)
        result = self._backend_executor.execute(
            command=command, arguments=args, check=True, location=tmp_dir
        )
        self._parse_output(tmp_dir)

        # self._remove_temporary(tmp_dir)
