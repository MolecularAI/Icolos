from typing import Dict, List
from icolos.core.containers.perturbation_map import Edge
from icolos.core.workflow_steps.pmx.base import StepPMXBase
from pydantic import BaseModel
from icolos.utils.enums.program_parameters import (
    GromacsEnum,
    StepPMXEnum,
)
from icolos.utils.enums.step_enums import StepGromacsEnum
from icolos.utils.execute_external.pmx import PMXExecutor
from icolos.utils.general.parallelization import SubtaskContainer
import os

_PSE = StepPMXEnum()
_SGE = StepGromacsEnum()
_GE = GromacsEnum()


class StepPMXPrepareSimulations(StepPMXBase, BaseModel):
    """
    Prepare the tpr file for either equilibration or production simulations
    """

    def __init__(self, **data):
        super().__init__(**data)

        self._initialize_backend(executor=PMXExecutor)

    def execute(self):

        edges = [e.get_edge_id() for e in self.get_edges()]
        self.execution.parallelization.max_length_sublists = 1
        self._subtask_container = SubtaskContainer(
            max_tries=self.execution.failure_policy.n_tries
        )
        self._subtask_container.load_data(edges)
        self._execute_pmx_step_parallel(
            run_func=self.prepare_simulation, step_id="pmx prepare_simulations"
        )

    def prepare_simulation(self, jobs: List[Edge], bLig=True, bProt=True):
        mdp_path = os.path.join(self.work_dir, "input/mdp")

        sim_type = self.settings.additional[_PSE.SIM_TYPE]

        for edge in jobs:

            ligTopPath = self._get_specific_path(
                workPath=self.work_dir, edge=edge, wp="water"
            )
            protTopPath = self._get_specific_path(
                workPath=self.work_dir, edge=edge, wp="protein"
            )

            for state in self.states:
                for r in range(1, self.get_perturbation_map().replicas + 1):

                    # ligand
                    if bLig == True:
                        wp = "water"
                        simpath = self._get_specific_path(
                            workPath=self.work_dir,
                            edge=edge,
                            wp=wp,
                            state=state,
                            r=r,
                            sim=sim_type,
                        )
                        empath = self._get_specific_path(
                            workPath=self.work_dir,
                            edge=edge,
                            wp=wp,
                            state=state,
                            r=r,
                            sim="em",
                        )
                        toppath = ligTopPath
                        self._prepare_single_tpr(
                            simpath, toppath, state, sim_type, empath
                        )

                    # protein
                    if bProt == True:
                        wp = "protein"
                        simpath = self._get_specific_path(
                            workPath=self.work_dir,
                            edge=edge,
                            wp=wp,
                            state=state,
                            r=r,
                            sim=sim_type,
                        )
                        empath = self._get_specific_path(
                            workPath=self.work_dir,
                            edge=edge,
                            wp=wp,
                            state=state,
                            r=r,
                            sim="em",
                        )
                        toppath = protTopPath
                        self._prepare_single_tpr(
                            simpath, toppath, state, sim_type, empath
                        )
