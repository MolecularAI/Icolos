from icolos.utils.enums.program_parameters import GromacsEnum, OpenBabelEnum
from icolos.utils.general.parallelization import Parallelizer, SubtaskContainer
from pydantic import BaseModel
from icolos.core.workflow_steps.fpocket.base import StepFpocketBase
from icolos.utils.enums.step_enums import StepCavExploreEnum, StepGromacsEnum
from icolos.utils.execute_external.execute import Executor
from icolos.utils.execute_external.openbabel import OpenBabelExecutor
from icolos.core.workflow_steps.step import _LE
from sklearn.cluster import DBSCAN
from collections import Counter
import numpy as np
import re
import os

_SFP = StepCavExploreEnum()
_OBE = OpenBabelEnum()
_SGE = StepGromacsEnum()


class StepMDpocket(StepFpocketBase, BaseModel):
    class Config:
        arbitrary_types_allowed = True

    obabel_executor: OpenBabelExecutor = None

    def __init__(self, **data):
        super().__init__(**data)

        self._initialize_backend(executor=Executor)
        self.obabel_executor = OpenBabelExecutor()

        # set max_length_sublists to 1
        self.execution.parallelization.max_length_sublists = 1

        pass

    def convert_gro_to_pdb(self, tmp_dir: str) -> None:
        structure = self.get_topol().structures[0]
        args = ["-igro", structure.get_file_name(), "-O", "confout.pdb"]
        self.obabel_executor.execute(
            command=_OBE.OBABEL, arguments=args, check=True, location=tmp_dir
        )
        os.remove(os.path.join(tmp_dir, structure.get_file_name()))

    def execute(self):
        pass
