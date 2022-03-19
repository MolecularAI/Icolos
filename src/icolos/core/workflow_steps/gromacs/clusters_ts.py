from icolos.core.workflow_steps.gromacs.base import StepGromacsBase
from icolos.utils.enums.step_enums import StepGromacsEnum
from typing import List
from pydantic import BaseModel
from icolos.utils.execute_external.execute import Executor
from icolos.core.workflow_steps.step import _LE
from icolos.utils.enums.program_parameters import GromacsEnum
import os
import sys

_SGE = StepGromacsEnum()
_GE = GromacsEnum()


class StepClusterTS(StepGromacsBase, BaseModel):
    """
    Generate time-resolved cluster plots from the output of gmx cluster, relies on MDplot R package
    """

    def __init__(self, **data):
        super().__init__(**data)

        self._initialize_backend(executor=Executor)

    def _construct_args(self, defaults: dict) -> List:
        """
        Custom method for argument construction, includes checks for required args in the config.
        """
        args = []
        for key, value in self.settings.arguments.parameters.items():
            args.append("".join([key, "=", value]))

        for value in self.settings.arguments.flags:
            args.append(value)
        for key, value in defaults.items():
            if key not in self.settings.arguments.parameters.keys():
                args.append("".join([key, "=", value]))

        # do some checks to make sure the required params have been passed
        for arg in [_SGE.CLUSTERS_NUMBER, _SGE.LENGTHS]:
            if arg not in self.settings.arguments.parameters.keys():
                self._logger.log(
                    f"Argument for parameter {arg} not found in provided argument. \
                        This must be specified!.  If this workflow has attached stdin, \
                            you can enter the value now...",
                    _LE.WARNING,
                )
                # instead of bailing out, take input from user if process has stdin connected
                if sys.stdin and sys.stdin.isatty():
                    value = input(f"Provide the parameter for option {arg}>>>")
                    args.append("".join([key, "=", value]))
                else:
                    self._logger.log(
                        f"No stdin stream detected, and cannot infer argument, step {self.step_id} may fail",
                        _LE.WARNING,
                    )
        return args

    def execute(self):
        """
        Visualise time-resolved gmx cluster results.
        Requires predceeding gmx_cluster step with clust-id.xvg file
        (ensure -clid flag is set, and xvg file is passed to this step)
        """

        tmp_dir = self._make_tmpdir()
        self.data.generic.write_out_all_files(tmp_dir)
        xvg_file = self.data.generic.get_argument_by_extension(ext="xvg")

        arguments = self._construct_args(
            defaults={
                "files": os.path.join(tmp_dir, xvg_file),
                "size": "1500,1500",
                "outformat": "png",
                "outfile": "clusters_ts.png",
                "timeUnit": "ns",
                "title": "CLUSTERS_timeseries",
            },
        )

        self._backend_executor.execute(
            command=_GE.CLUSTER_TS, arguments=arguments, location=tmp_dir, check=True
        )

        self._parse_output(tmp_dir)
        self._remove_temporary(tmp_dir)
