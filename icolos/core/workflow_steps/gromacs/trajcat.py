from icolos.utils.enums.program_parameters import GromacsEnum
from icolos.core.workflow_steps.gromacs.base import StepGromacsBase
from icolos.utils.execute_external.gromacs import GromacsExecutor
from pydantic import BaseModel
from icolos.core.workflow_steps.step import _LE
import os

_GE = GromacsEnum()


class StepGMXTrjcat(StepGromacsBase, BaseModel):
    """
    Concatenates multiple trajectories, useful for subsequent rmsd/cluster calculations
    """

    def __init__(self, **data):
        super().__init__(**data)

        self._initialize_backend(executor=GromacsExecutor)
        self._check_backend_availability()

    def execute(self):

        tmp_dir = self._make_tmpdir()

        # write the trajectories to the tmpdir, writing to separate file names, then glob the xtc files

        for idx, file in enumerate(self.data.generic.get_files_by_extension(ext="xtc")):
            file.write(path=os.path.join(tmp_dir, f"traj_{idx}.xtc"), join=False)

        flag_dict = {
            "-f": "*.xtc",
            "-o": "trjcat_out.xtc",
            "-cat": "",  # need this to paste the trajectories back to back
        }

        arguments = self._parse_arguments(flag_dict=flag_dict)
        result = self._backend_executor.execute(
            command=_GE.TRJCAT, arguments=arguments, location=tmp_dir, check=True
        )
        for line in result.stdout.split("\n"):
            self._logger_blank.log(line, _LE.DEBUG)

        rm_files = [
            f for f in os.listdir(tmp_dir) if f.endswith("xtc") and "trjcat" not in f
        ]
        for f in rm_files:
            os.remove(os.path.join(tmp_dir, f))
        self._parse_output(tmp_dir)
        self._remove_temporary(tmp_dir)
