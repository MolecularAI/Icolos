# implement pdbfixer as FOSS alternative to proteinprep
from icolos.utils.enums.step_enums import StepPdbFixerEnum
from icolos.core.workflow_steps.step import StepBase
from icolos.utils.enums.program_parameters import PdbFixerEnum
from icolos.utils.execute_external.execute import Executor
from pydantic import BaseModel
from pdbfixer.pdbfixer import PDBFixer
import os


_SFE = StepPdbFixerEnum()
_FE = PdbFixerEnum()


class StepPdbFixer(StepBase, BaseModel):
    def __init__(self, **data):
        super().__init__(**data)
        self._initialize_backend(executor=Executor)

    def _parse_arguments(self):
        default_flags = [
            "--replace-nonstandard",
            "--add-residues",
        ]
        default_params = {
            "--ph": "7.0",
            "--add-atoms": "all",
            "--keep-heterogens": "all",
        }
        arguments = []
        for arg in self.settings.arguments.flags:
            arguments.append(arg)
        for key, value in self.settings.arguments.parameters.items():
            formatted_arg = f"{key}={value}"
            arguments.append(formatted_arg)
        for key in default_flags:
            if key not in self.settings.arguments.flags:
                arguments.append(key)
        for key, value in default_params.items():
            if key not in self.settings.arguments.parameters.keys():
                formatted_arg = f"{key}={value}"
                arguments.append(formatted_arg)
        return arguments

    def execute(self):

        tmp_dir = self._make_tmpdir()

        self.data.generic.write_out_all_files(tmp_dir)
        pdb_files = self.data.generic.get_file_names_by_extension("pdb")

        arguments = self._parse_arguments()

        for file in pdb_files:
            path = os.path.join(tmp_dir, file)
            arguments.extend(["--output", path])
            arguments = [path] + arguments

            self._backend_executor.execute(
                command=_FE.FIXER, arguments=arguments, location=tmp_dir, check=True
            )

            #
        self._parse_output(tmp_dir)

        self._remove_temporary(tmp_dir)
