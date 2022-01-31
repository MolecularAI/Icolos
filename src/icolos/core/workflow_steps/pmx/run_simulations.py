from typing import Dict, List
from icolos.core.workflow_steps.pmx.base import StepPMXBase
from pydantic import BaseModel
import numpy as np
from icolos.utils.enums.program_parameters import PMXAtomMappingEnum, PMXEnum
from icolos.utils.execute_external.pmx import PMXExecutor
from icolos.utils.general.parallelization import SubtaskContainer

_PE = PMXEnum()
_PAE = PMXAtomMappingEnum()


class StepPMXRunSimulations(StepPMXBase, BaseModel):
    """
    Calls pmx run_simulations entrypoint, handles parallel execution across multiple GPUs
    """

    def __init__(self, **data):
        super().__init__(**data)

        self._initialize_backend(executor=PMXExecutor)

    def execute(self):

        edges = self.get_edges()
        # run everything through in one batch, with multiple edges per call
        self.execution.parallelization.max_length_sublists = int(
            np.ceil(len(edges) / self._get_number_cores())
        )
        self._subtask_container = SubtaskContainer(
            max_tries=self.execution.failure_policy.n_tries
        )
        self._subtask_container.load_data(edges)
        self._execute_pmx_step_parallel(
            run_func=self._execute_command, step_id="pmx_run_simulations"
        )

    def _execute_command(self, edges: List, q: Dict):
        """
        Execute the simulations for one edge, calling the run_sims pmx entrypoint
        """
        args = {
            "-edges": '"' + " ".join([e.get_edge_id() for e in edges]) + '"',
            "-workPath": self.work_dir,
            "-sim_type": self.settings.additional["sim_type"],
            "-replicas": self.get_workflow_object().workflow_data.perturbation_map.replicas,
        }
        for key, value in self.settings.arguments.parameters:
            args[key] = value

        result = self._backend_executor.execute(
            command=_PE.RUN_SIMULATIONS,
            arguments=self.get_arguments(defaults=args),
            check=True,
            location=self.work_dir,
        )

        q[edges[0].get_edge_id()] = result.returncode
