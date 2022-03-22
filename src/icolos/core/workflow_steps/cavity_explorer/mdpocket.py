from icolos.utils.general.parallelization import Parallelizer, SubtaskContainer
from pydantic import BaseModel
from icolos.core.workflow_steps.cavity_explorer.base import StepCavityExplorerBase
from icolos.utils.enums.step_enums import StepCavExploreEnum
from icolos.utils.execute_external.execute import Executor
from icolos.core.workflow_steps.step import _LE
from sklearn.cluster import DBSCAN
from collections import Counter
import numpy as np
import re
import os

_SFP = StepCavExploreEnum()


class StepMDpocket(StepCavityExplorerBase, BaseModel):
    def __init__(self, **data):
        super().__init__(**data)

        # self._initialize_backend(executor=MPIExecutor)
        self._initialize_backend(executor=Executor)

        # set max_length_sublists to 1
        self.execution.parallelization.max_length_sublists = 1

        pass

    def execute(self):
        pass
