from icolos.utils.enums.logging_enums import LoggingConfigEnum
from typing import Dict, List
from icolos.core.workflow_steps.pmx.base import StepPMXBase
from pydantic import BaseModel
from icolos.utils.enums.program_parameters import PMXAtomMappingEnum, PMXEnum
from icolos.utils.execute_external.pmx import PMXExecutor
from icolos.utils.general.parallelization import SubtaskContainer

_PE = PMXEnum()
_PAE = PMXAtomMappingEnum()
_LE = LoggingConfigEnum()


class StepPMXPrepareTransitions(StepPMXBase, BaseModel):
    """
    Executes the pmx prepare_transitions.py entrypoint
    """

    def __init__(self, **data):
        super().__init__(**data)

        self._initialize_backend(executor=PMXExecutor)

    def execute(self):
        edges = self.get_edges()

        self.execution.parallelization.max_length_sublists = 1
        self._subtask_container = SubtaskContainer(
            max_tries=self.execution.failure_policy.n_tries
        )
        self._subtask_container.load_data(edges)
        self._execute_pmx_step_parallel(
            run_func=self._execute_command, step_id="pmx_prepare_transitions"
        )

    def _execute_command(self, edges: List, q: Dict):
        args = {
            "-edges": '"' + " ".join([e.get_edge_id() for e in edges]) + '"',
            "-workPath": self.work_dir,
            "-replicas": self.get_workflow_object().workflow_data.perturbation_map.replicas,
        }
        result = self._backend_executor.execute(
            command=_PE.PREPARE_TRANSITIONS,
            arguments=self.get_arguments(defaults=args),
            check=True,
            location=self.work_dir,
        )
        q[edges[0].get_edge_id()] = result.returncode
