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


class StepGMXGenion(StepGromacsBase, BaseModel):
    """
    Wrapper for gmx genion
    """

    def __init__(self, **data):
        super().__init__(**data)

        self._initialize_backend(executor=GromacsExecutor)
        self._check_backend_availability()

    def execute(self):
        tmp_dir = self._make_tmpdir()
        topol = self.get_topol()
        self.write_input_files(tmp_dir, topol=topol)
        arguments = self._parse_arguments(
            {
                "-o": _SGE.STD_STRUCTURE,
                "-p": _SGE.STD_TOPOL,
                "-s": _SGE.STD_TPR,
            }
        )
        result = self._backend_executor.execute(
            command=_GE.GENION,
            arguments=arguments,
            location=tmp_dir,
            pipe_input=self.construct_pipe_arguments(
                tmp_dir, self.settings.additional[_SBE.PIPE_INPUT]
            ),
        )
        for line in result.stdout.split("\n"):
            self._logger_blank.log(line, _LE.DEBUG)
        self._logger.log(
            f"Completed execution for {self.step_id} successfully", _LE.INFO
        )
        # this is the last structural change to the topology in a regular gromacs setup,
        # update the index groups here
        make_ndx_args = ["-f", _SGE.STD_STRUCTURE, "-o", _SGE.STD_INDEX]
        index_files = [f for f in os.listdir(tmp_dir) if f.endswith(".ndx")]
        # remove any existing index files
        for f in index_files:
            self._remove_temporary(os.path.join(tmp_dir, f))
        # generate new index file
        result = self._backend_executor.execute(
            command=_GE.MAKE_NDX,
            arguments=make_ndx_args,
            location=tmp_dir,
            check=True,
            pipe_input='echo -e "1 | 12 \nq"',
        )
        for line in result.stdout.split("\n"):
            self._logger_blank.log(line, _LE.DEBUG)

        self._logger.log('Added index group to "index.ndx"', _LE.DEBUG)
        self._parse_output(tmp_dir)
        topol.set_structure(tmp_dir)
        topol.parse(tmp_dir)
        self._remove_temporary(tmp_dir)
