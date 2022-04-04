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


class StepGMXCluster(StepGromacsBase, BaseModel):
    """
    Execute gmx cluster on a trajectory
    """

    def __init__(self, **data):
        super().__init__(**data)

        self._initialize_backend(executor=GromacsExecutor)
        self._check_backend_availability()

    def execute(self):
        tmp_dir = self._make_tmpdir()
        topol = self.get_topol()
        self.write_input_files(tmp_dir, topol=topol)

        # give the option to run a make_ndx step preceding clustering to facilitate clustering on custom groups
        if _SGE.INDEX_FLAG in self.settings.arguments.parameters.keys():
            assert (
                _SGE.STD_INDEX in os.listdir(tmp_dir)
                or self.settings.additional[_SGE.MAKE_NDX_COMMAND] is not None
            )
            if _SGE.STD_INDEX not in os.listdir(tmp_dir):
                try:
                    ndx_arguments = [
                        "-f",
                        _SGE.STD_STRUCTURE,
                        "-o",
                        _SGE.STD_INDEX,
                    ]
                    result = self._backend_executor.execute(
                        command=_GE.MAKE_NDX,
                        arguments=ndx_arguments,
                        location=tmp_dir,
                        check=True,
                        pipe_input=self.construct_pipe_arguments(
                            tmp_dir, self.settings.additional[_SGE.MAKE_NDX_COMMAND]
                        ),
                    )

                except KeyError:
                    raise KeyError(
                        "If the index flag was specified, you must provide the ndx command in additional "
                        "settings"
                    )

        flag_dict = {
            "-s": _SGE.STD_TPR,
            "-f": _SGE.STD_XTC,
            "-cl": "clusters.pdb",
            "-clid": "cluster_id.xvg",
        }
        arguments = self._parse_arguments(flag_dict=flag_dict)

        result = self._backend_executor.execute(
            command=_GE.CLUSTER,
            arguments=arguments,
            location=tmp_dir,
            check=True,
            pipe_input=self.construct_pipe_arguments(
                tmp_dir, self.settings.additional[_SBE.PIPE_INPUT]
            ),
        )

        for line in result.stdout.split("\n"):
            self._logger_blank.log(line, _LE.DEBUG)
        self._parse_output(tmp_dir)
        self._remove_temporary(tmp_dir)
