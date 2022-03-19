from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum
from icolos.utils.enums.program_parameters import GromacsEnum
from icolos.core.workflow_steps.gromacs.base import StepGromacsBase
from icolos.utils.execute_external.gromacs import GromacsExecutor
from pydantic import BaseModel
from icolos.core.workflow_steps.step import _LE
import os

_GE = GromacsEnum()
_SGE = StepGromacsEnum()
_SBE = StepBaseEnum


class StepGMXTrjconv(StepGromacsBase, BaseModel):
    """
    Postprocessing step for gromacs trajectories
    Mostly used for removing pbc, fitting trajectory etc.
    """

    def __init__(self, **data):
        super().__init__(**data)

        self._initialize_backend(executor=GromacsExecutor)
        self._check_backend_availability()

    def execute(self):

        tmp_dir = self._make_tmpdir()
        topol = self.get_topol()
        topol.write_ndx(tmp_dir)
        for i in range(len(topol.trajectories.keys())):
            topol.write_tpr(tmp_dir, index=i)
            topol.write_trajectory(tmp_dir, index=i)
            fitted_traj = f"traj_{i}.xtc"
            flag_dict = {
                "-s": _SGE.STD_TPR,
                "-f": _SGE.STD_XTC,
                "-o": fitted_traj,
                "-n": _SGE.STD_INDEX,
            }
            arguments = self._parse_arguments(flag_dict=flag_dict)
            result = self._backend_executor.execute(
                command=_GE.TRJCONV,
                arguments=arguments,
                location=tmp_dir,
                check=True,
                pipe_input=self.construct_pipe_arguments(
                    tmp_dir, self.settings.additional[_SBE.PIPE_INPUT]
                ),
            )
            # ensure the old trajectories are not parsed for writeout
            os.remove(os.path.join(tmp_dir, _SGE.STD_XTC))
            for line in result.stdout.split("\n"):
                self._logger_blank.log(line, _LE.DEBUG)
            topol.set_trajectory(path=tmp_dir, file=fitted_traj, index=i)
            topol.set_structure(path=tmp_dir, index=i)
        self._parse_output(tmp_dir)
        self._remove_temporary(tmp_dir)
