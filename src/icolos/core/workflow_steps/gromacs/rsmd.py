from icolos.utils.enums.program_parameters import GromacsEnum
from icolos.core.workflow_steps.gromacs.base import StepGromacsBase
from icolos.utils.execute_external.gromacs import GromacsExecutor
from pydantic import BaseModel
from icolos.core.workflow_steps.step import _LE
import os

_GE = GromacsEnum()


class StepGMXrmsd(StepGromacsBase, BaseModel):
    """
    Run gromacs rmsd calculation on trajectory
    """

    def __init__(self, **data):
        super().__init__(**data)

        self._initialize_backend(executor=GromacsExecutor)
        self._check_backend_availability()

    def execute(self):

        tmp_dir = self._make_tmpdir()

        # we're going to get a trajectory from a trjconv step, and a single structure to
        # compare against. rmsd is computed for every step of the trj file

        # write out generic files
        self.write_input_files(tmp_dir)

        # conformer coming from a Compound object
        conf = self._unroll_compounds(self.data.compounds)

        conf = conf[0]
        conf.write(os.path.join(tmp_dir, "reference.sdf"), format_="pdb")

        flag_dict = {
            "-s": "reference.pdb",
            "-f": self.data.generic.get_argument_by_extension("xtc"),
            "-fit": "rot+trans",
        }

        arguments = self._parse_arguments(flag_dict=flag_dict, args=["-w"])
        result = self._backend_executor.execute(
            command=_GE.RMS,
            arguments=arguments,
            location=tmp_dir,
            check=True,
            pipe_input='echo -e "2\n2\n"',
        )
        for line in result.stdout.split("\n"):
            self._logger_blank.log(line, _LE.DEBUG)

        self._parse_output(tmp_dir)
        self._remove_temporary(tmp_dir)
