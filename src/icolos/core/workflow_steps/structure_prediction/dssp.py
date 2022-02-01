from typing import List
from icolos.core.workflow_steps.step import StepBase
from icolos.utils.execute_external.execute import Executor
from pydantic import BaseModel
from icolos.utils.enums.step_enums import StepDSSPEnum
from icolos.utils.enums.program_parameters import DSSPEnum
import os


_SDE = StepDSSPEnum()
_DE = DSSPEnum()


class StepDSSP(StepBase, BaseModel):
    def __init__(self, **data):
        super().__init__(**data)
        self._initialize_backend(executor=Executor)

    def _construct_arguments(self, tmp_dir: str, file: str) -> List:
        args = []
        for flag in self.settings.arguments.flags:
            args.append(flag)
        for key, value in self.settings.arguments.parameters.items():
            args.append(key)
            args.append(value)

        # set the input and output files
        args.append(file)
        output = f"dssp_output_{file.split('.')[0]}.txt"
        args.append(output)
        return args

    def _parse_output(self, tmp_dir: str) -> None:
        for file in [f for f in os.listdir(tmp_dir) if f.endswith("txt")]:
            with open(os.path.join(tmp_dir, file), "r") as f:
                self._add_data_to_generic(file, f.read())

    def execute(self):
        """
        Executes dssp on a set of input structures
        """

        tmp_dir = self._make_tmpdir()
        print(tmp_dir)
        self.data.generic.write_out_all_files(tmp_dir)

        file_list = self.data.generic.get_file_names_by_extension(ext="pdb")

        for file in file_list:
            arguments = self._construct_arguments(tmp_dir, file)
            self._backend_executor.execute(
                command=_DE.MKDSSP, arguments=arguments, check=True, location=tmp_dir
            )

        self._parse_output(tmp_dir)
        self._remove_temporary(tmp_dir)
