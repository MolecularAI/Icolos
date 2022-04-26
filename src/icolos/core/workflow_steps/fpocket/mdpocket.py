from icolos.utils.general.parallelization import Parallelizer, SubtaskContainer
from pydantic import BaseModel
from icolos.core.workflow_steps.fpocket.base import StepFpocketBase
from icolos.utils.enums.step_enums import StepCavExploreEnum
from icolos.utils.execute_external.execute import Executor
from icolos.core.workflow_steps.step import _LE
from sklearn.cluster import DBSCAN
from collections import Counter
import numpy as np
import re
import os

_SFP = StepCavExploreEnum()


class StepMDpocket(StepFpocketBase, BaseModel):
    def __init__(self, **data):
        super().__init__(**data)

        # self._initialize_backend(executor=MPIExecutor)
        self._initialize_backend(executor=Executor)

        # set max_length_sublists to 1
        self.execution.parallelization.max_length_sublists = 1

    def _create_density_grid_file(self, tmp_dir: str, iso_value: float):
        """creates a density grid from the .dx-file into a .pdb-file, heavily influenced by extractISOPdb.py provided
        by fpocket"""
        pass

    def _cluster_pockets(self, tmp_dir, eps, min_samples, threshold, iso_value):
        """
        Clusters points from the initial MDpocket density grid, at a certain iso value
        """
        pass

    def _save_pocket_files(self, tmp_dir, data, labels):
        """saves the individual pockets as individual pdbs to be used with mdpocket"""
        pass

    def _run_mdpocket_selected_pocket(self, tmp_dir):
        """runs the second mdpocket command for fpocket3"""
        pass

    def _prepare_batch_inputs(self, batch, tmp_dir):
        tmp_dirs = []
        args = []
        for next_subtask_list in batch:
            tmp_dirs.append(tmp_dir)
            for (
                subtask
            ) in (
                next_subtask_list
            ):  # enforced only one task per subtask, otherwise it makes no sense
                args.append(subtask.data)  # append the arguments list
        return tmp_dirs, args

    def _execute_mdpocket(self, tmp_dir, arguments):

        self._backend_executor.execute(
            command=_SFP.MDPOCKET_COMMAND,
            arguments=arguments,
            location=tmp_dir,
            check=True,
        )

    def execute(self):
        pass
        