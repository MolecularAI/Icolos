from typing import Dict, List
from icolos.core.workflow_steps.pmx.base import StepPMXBase
from pydantic import BaseModel
import os
from icolos.utils.enums.program_parameters import PMXAtomMappingEnum, PMXEnum
from icolos.core.workflow_steps.step import _LE
from icolos.utils.execute_external.pmx import PMXExecutor
from icolos.utils.general.parallelization import SubtaskContainer

_PE = PMXEnum()
_PAE = PMXAtomMappingEnum()


class StepPMXBoxWaterIons(StepPMXBase, BaseModel):
    """
    Take the prepard structure files and prepare the system,
    runs editconf, solvate, genion and grompp for each system
    to be simulated
    """

    # Note all paths are relative to the workdir
    def __init__(self, **data):
        super().__init__(**data)

        self._initialize_backend(executor=PMXExecutor)

    def execute(self):
        # run the wrapper script in pmx to prepare the systems

        edges = self.get_edges()

        self.execution.parallelization.max_length_sublists = 1
        self._subtask_container = SubtaskContainer(
            max_tries=self.execution.failure_policy.n_tries
        )
        self._subtask_container.load_data(edges)
        self._execute_pmx_step_parallel(
            run_func=self._execute_command, step_id="BoxWaterIons"
        )

    def _execute_command(self, edges: List, q: Dict):

        arguments = {
            "-edges": '"' + " ".join([e.get_edge_id() for e in edges]) + '"',
            "-ligandPath": os.path.join(self.work_dir, _PAE.LIGAND_DIR),
            "-workPath": self.work_dir,
        }

        result = self._backend_executor.execute(
            command=_PE.BOX_WATER_IONS,
            arguments=self.get_arguments(defaults=arguments),
            check=True,
            location=self.work_dir,
        )

        self._logger.log("End of BoxWaterIons output", _LE.DEBUG)
        # collect returncodes from subprocess
        q[edges[0].get_edge_id()] = result.returncode
