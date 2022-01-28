from typing import Dict, List
from icolos.core.workflow_steps.pmx.base import StepPMXBase
from pydantic import BaseModel
import os
from icolos.utils.enums.program_parameters import PMXAtomMappingEnum, PMXEnum
from icolos.utils.execute_external.pmx import PMXExecutor
from icolos.utils.general.parallelization import SubtaskContainer

_PE = PMXEnum()
_PAE = PMXAtomMappingEnum()


class StepPMXAssembleSystems(StepPMXBase, BaseModel):
    """
    Executes the assemble_systems.py script, edges are parallelized over available cores

    """

    def __init__(self, **data):
        super().__init__(**data)

        self._initialize_backend(executor=PMXExecutor)

    def execute(self):
        assert self.work_dir is not None and os.path.isdir(self.work_dir)

        # get edges from the perturbation map attached to the step
        edges = self.get_edges()

        # enforce one edge per task list (results in multiple batches for large maps)
        self.execution.parallelization.max_length_sublists = 1
        self._subtask_container = SubtaskContainer(
            max_tries=self.execution.failure_policy.n_tries
        )
        self._subtask_container.load_data(edges)
        self._execute_pmx_step_parallel(
            run_func=self._execute_command, step_id="pmx assemble_systems"
        )

    def _execute_command(self, jobs: List):

        args = {
            "-edges": '"' + " ".join([e.get_edge_id() for e in jobs]) + '"',
            "-ligand_path": os.path.join(self.work_dir, _PAE.LIGAND_DIR),
            "-workPath": self.work_dir,
        }
        self._backend_executor.execute(
            command=_PE.ASSEMBLE_SYSTEMS,
            arguments=self.get_arguments(defaults=args),
            check=True,
            location=self.work_dir,
        )
