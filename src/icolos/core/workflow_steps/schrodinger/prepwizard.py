from icolos.utils.enums.step_enums import StepGromacsEnum, StepPrepwizEnum
from icolos.utils.enums.program_parameters import SchrodingerExecutablesEnum
from icolos.core.workflow_steps.schrodinger.base import StepSchrodingerBase
from icolos.utils.execute_external.schrodinger import SchrodingerExecutor
from icolos.core.containers.generic import GenericData
from pydantic import BaseModel
from copy import deepcopy
import os

_SEE = SchrodingerExecutablesEnum()
_SGE = StepGromacsEnum()
_SPE = StepPrepwizEnum()


class StepPrepwizard(StepSchrodingerBase, BaseModel):
    """
    Interface to Schrodinger's PrepWizard program for protein prep
    """

    def __init__(self, **data):
        super().__init__(**data)
        self._initialize_backend(executor=SchrodingerExecutor)
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
        input_file = self.data.generic.get_file_names_by_extension("pdb")[0]
        output_file = input_file  # write to the same file name to keep things tidy
        arguments.append(input_file)
        arguments.append(output_file)
        return arguments

    def _parse_output(self, tmp_dir: str):
        output_pdb = os.path.join(
            tmp_dir, self.data.generic.get_file_names_by_extension("pdb")[0]
        )
        with open(output_pdb, "r") as f:
            data = f.read()
        self.data.generic.clear_file_dict()
        output_file = GenericData(file_name=_SGE.COMPLEX_PDB, file_data=data)
        self.data.generic.add_file(output_file)

    def _remove_ligand(self, tmp_dir):
        remove_res = self.settings.additional[_SPE.REMOVE_RES]
        pdb_file = self.data.generic.get_argument_by_extension("pdb")
        cleaned_pdb_lines = []
        # handle ligand removal mode: strip ligands, leave cofactors
        if remove_res != _SPE.LIGANDS and not isinstance(remove_res, list):
            remove_res = list(remove_res)

        with open(os.path.join(tmp_dir, pdb_file), "r") as f:
            if remove_res == _SPE.LIGANDS:
                # automatically remove ligands, keep cofactors that are specified in the enum.
                for line in f.readlines():
                    if (
                        line is not None
                        and len(line.split()) > 3
                        and (
                            line.split()[0] == "ATOM"
                            or any(l in line for l in _SPE.COFACTOR_IDS)
                        )
                    ):
                        cleaned_pdb_lines.append(line)
            else:
                for line in f.readlines():
                    if not any(l in line for l in remove_res):
                        cleaned_pdb_lines.append(line)

        with open(os.path.join(tmp_dir, pdb_file), "w") as f:
            f.writelines(cleaned_pdb_lines)

    def execute(self):
        tmp_dir = self._make_tmpdir()
        args = self._parse_args()
        self.data.generic.write_out_all_files(tmp_dir)
        if (
            _SPE.REMOVE_RES in self.settings.additional.keys()
            and self.settings.additional[_SPE.REMOVE_RES] is not None
        ):
            self._remove_ligand(tmp_dir)
        self._backend_executor.execute(
            command=_SEE.PREPWIZARD, arguments=args, check=True, location=tmp_dir
        )

        self._parse_output(tmp_dir)
