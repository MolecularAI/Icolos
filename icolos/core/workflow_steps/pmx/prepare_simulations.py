from typing import Dict, List
from icolos.core.workflow_steps.pmx.base import StepPMXBase
from pydantic import BaseModel
from icolos.utils.enums.program_parameters import PMXAtomMappingEnum, PMXEnum
from icolos.utils.execute_external.pmx import PMXExecutor
from icolos.utils.general.parallelization import SubtaskContainer

_PE = PMXEnum()
_PAE = PMXAtomMappingEnum()


class StepPMXPrepareSimulations(StepPMXBase, BaseModel):
    """
    Prepare the tpr file for either equilibration or production simulations

    Calls pmx util entrypoint prepare_simulations.py with
    list of edges and the workdir path
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
            run_func=self._execute_command, step_id="pmx_prepare_sims"
        )

    def _execute_command(self, edges: List, q: Dict):
        arguments = {
            "-edges": '"' + " ".join([e.get_edge_id() for e in edges]) + '"',
            "-workPath": self.work_dir,
            "-sim_type": self.settings.additional["sim_type"],
            "-replicas": self.get_workflow_object().workflow_data.perturbation_map.replicas,
        }
        result = self._backend_executor.execute(
            command=_PE.PREPARE_SIMULATIONS,
            arguments=self.get_arguments(defaults=arguments),
            check=True,
            location=self.work_dir,
        )

        q[edges[0].get_edge_id()] = result.returncode
