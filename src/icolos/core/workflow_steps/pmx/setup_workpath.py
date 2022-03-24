import os
from typing import List
from pydantic import BaseModel
from icolos.core.containers.perturbation_map import Node
from icolos.core.workflow_steps.pmx.base import StepPMXBase
from icolos.utils.enums.program_parameters import GromacsEnum
from icolos.utils.enums.step_enums import StepGromacsEnum
from icolos.utils.execute_external.execute import Executor
from icolos.utils.execute_external.gromacs import GromacsExecutor
from icolos.utils.general.parallelization import SubtaskContainer

_GE = GromacsEnum()
_SGE = StepGromacsEnum()

# These classes are based on the work of Vytautas Gapsys et al: https://github.com/deGrootLab/pmx/
class StepPMXSetup(StepPMXBase, BaseModel):
    """
    Create the directory tree structure.
    Requires the pmx workflow to be executing using the single_dir running mode
    Operates on the perturbation map object, runs acpype
    on the written structures to produce the amber-compatible itp files
    """

    _gromacs_executor: GromacsExecutor = None
    _antechamber_executor: Executor = None

    def __init__(self, **data):
        super().__init__(**data)
        self._gromacs_executor = GromacsExecutor(
            prefix_execution=self.execution.prefix_execution
        )

    def execute(self):
        # sets the number of replicas to be used throughput the pmx run
        replicas = (
            self.settings.additional["replicas"]
            if "replicas" in self.settings.additional.keys()
            else 3
        )
        assert self.work_dir is not None and os.path.isdir(self.work_dir)
        self._construct_perturbation_map(self.work_dir, replicas)
        # create the directory structure for subsequent calculations
        edges = self.get_edges()
        nodes = self.get_nodes()

        # create the input directory to sit at the top level of the workdir, contains ligands,
        # mdp and protein topology files
        os.makedirs(os.path.join(self.work_dir, "input"), exist_ok=True)
        for folder in ["ligands", "mdp", "protein"]:
            os.makedirs(os.path.join(self.work_dir, "input", folder), exist_ok=True)

        # handle protein parametrization with pdb2gmx
        protein = (
            self.get_workflow_object().workflow_data.perturbation_map.get_protein()
        )
        protein.write(os.path.join(self.work_dir, "input/protein"))

        self._parametrise_protein(protein=protein.get_file_name(), path="input/protein")

        # remove the backup file
        old_protein = [
            f
            for f in os.listdir(os.path.join(self.work_dir, "input/protein"))
            if f.endswith("#")
        ]
        # only want the parametrised processed pdb file in there
        old_protein.append(protein.get_file_name())
        for f in old_protein:
            os.remove(os.path.join(self.work_dir, "input/protein", f))

        self._clean_protein()

        mdp_dir = self.data.generic.get_argument_by_extension(
            ext="mdp", rtn_file_object=True
        )
        mdp_dir.write(os.path.join(self.work_dir, "input/mdp"))

        # parallelize the antechamber call across the pool of nodes

        self.execution.parallelization.max_length_sublists = 1
        self._subtask_container = SubtaskContainer(
            max_tries=self.execution.failure_policy.n_tries
        )
        self._subtask_container.load_data(nodes)
        self._execute_pmx_step_parallel(
            run_func=self._parametrise_nodes,
            step_id="pmx_setup",
            result_checker=self._check_results,
        )

        # create the output folder structure
        for edge in edges:
            edgepath = os.path.join(
                self.work_dir,
                str(f"{edge.node_from.get_node_hash()}_{edge.node_to.get_node_hash()}"),
            )
            hybridTopFolder = f"{edgepath}/hybridStrTop"
            os.makedirs(hybridTopFolder, exist_ok=True)

            # water/protein
            for wp in self.therm_cycle_branches:
                wppath = f"{edgepath}/{wp}"
                os.makedirs(wppath, exist_ok=True)

                # stateA/stateB
                for state in self.states:
                    statepath = f"{wppath}/{state}"
                    os.makedirs(statepath, exist_ok=True)

                    # run1/run2/run3
                    for r in range(1, replicas + 1):
                        runpath = f"{statepath}/run{r}"
                        os.makedirs(runpath, exist_ok=True)

                        # em/eq_posre/eq/transitions
                        for sim in self.sim_types:
                            simpath = f"{runpath}/{sim}".format(runpath, sim)
                            os.makedirs(simpath, exist_ok=True)

    def _check_results(self, batch: List[List[Node]]) -> List[List[bool]]:
        output_files = ["ffMOL.itp", "MOL.itp"]
        results = []
        for subjob in batch:
            subjob_results = []
            for job in subjob:
                subjob_results.append(
                    all(
                        [
                            os.path.isfile(
                                os.path.join(
                                    self.work_dir,
                                    "input/ligands",
                                    job.get_node_hash(),
                                    f,
                                )
                            )
                            for f in output_files
                        ]
                    )
                )
            results.append(subjob_results)
        return results
