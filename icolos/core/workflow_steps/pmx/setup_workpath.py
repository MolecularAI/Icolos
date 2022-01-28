from icolos.core.containers.perturbation_map import Node
import os
from typing import Dict
from pydantic import BaseModel
from icolos.core.workflow_steps.pmx.base import StepPMXBase
from icolos.utils.enums.program_parameters import GromacsEnum
from icolos.utils.enums.step_enums import StepGromacsEnum
from icolos.utils.execute_external.execute import Executor
from icolos.utils.execute_external.gromacs import GromacsExecutor
from icolos.utils.general.parallelization import SubtaskContainer

_GE = GromacsEnum()
_SGE = StepGromacsEnum()


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
        self._gromacs_executor = GromacsExecutor(prefix_execution=_SGE.GROMACS_LOAD)

    def _separate_atomtypes(self, lig_path: str) -> None:
        with open(os.path.join(lig_path, "MOL.itp"), "r") as f:
            itp_lines = f.readlines()

        start_idx = self._get_line_idx(itp_lines, _GE.ATOMTYPES)
        stop_index = self._get_line_idx(itp_lines, _GE.MOLECULETYPES)

        atomtype_lines = itp_lines[start_idx:stop_index]
        cleaned_itp_lines = itp_lines[stop_index:]
        with open(os.path.join(lig_path, "MOL.itp"), "w") as f:
            f.writelines(cleaned_itp_lines)

        # process the atomtype lines to remove the bondtype
        # col causes gmx to complain
        cleaned_atomtype_lines = []
        for line in atomtype_lines:
            parts = line.split()
            if len(parts) > 5:
                cleaned_parts = [parts[0]] + parts[2:] + ["\n"]
                cleaned_atomtype_lines.append(" ".join(cleaned_parts))
        with open(os.path.join(lig_path, "ffMOL.itp"), "w") as f:
            f.writelines(cleaned_atomtype_lines)

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

        existing_itp_files = [
            f
            for f in os.listdir(os.path.join(self.work_dir, "input/protein"))
            if f.endswith("itp") and f.startswith("Protein")
        ]
        if (
            not existing_itp_files
        ):  # no protein itp files, we have a single chain that needs extacting from the top file
            with open(os.path.join(self.work_dir, "input/protein/topol.top"), "r") as f:
                top_lines = f.readlines()

            moltype_line = self._get_line_idx(top_lines, _GE.MOLECULETYPES)

            end_itp_line = self._get_line_idx(top_lines, "; Include water topology")

            moltype = top_lines[moltype_line + 2].split()[0]
            cleaned_top = (
                top_lines[:moltype_line]
                + [f'#include "topol_{moltype}.itp']
                + top_lines[end_itp_line:]
            )

            itp_lines = top_lines[moltype_line:end_itp_line]

            with open(os.path.join(self.work_dir, "input/protein/topol.top"), "w") as f:
                f.writelines(cleaned_top)

            with open(
                os.path.join(self.work_dir, f"input/protein/topol_{moltype}.itp"), "w"
            ) as f:
                f.writelines(itp_lines)

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
            run_func=self._parametrise_nodes, step_id="pmx_setup"
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
            for wp in ["water", "protein"]:
                wppath = f"{edgepath}/{wp}"
                os.makedirs(wppath, exist_ok=True)

                # stateA/stateB
                for state in ["stateA", "stateB"]:
                    statepath = f"{wppath}/{state}"
                    os.makedirs(statepath, exist_ok=True)

                    # run1/run2/run3
                    for r in range(1, replicas + 1):
                        runpath = f"{statepath}/run{r}"
                        os.makedirs(runpath, exist_ok=True)

                        # em/eq_posre/eq/transitions
                        for sim in ["em", "eq", "transitions"]:
                            simpath = f"{runpath}/{sim}".format(runpath, sim)
                            os.makedirs(simpath, exist_ok=True)

    # TODO: sort out nomenclature here
    def _parametrise_nodes(self, edges: Node, q: Dict):
        # because we use the base-class infrastructure to parallelize, arg names are awkward
        # in this case, we parallize over nodes, not edges!
        if isinstance(edges, list):
            node = edges[0]
        else:
            node = edges
        lig_path = os.path.join(self.work_dir, "input", "ligands", node.get_node_hash())
        os.makedirs(lig_path, exist_ok=True)
        node.conformer.write(os.path.join(lig_path, "MOL.sdf"), format_="pdb")

        # clean the written pdb, remove anything except hetatm/atom lines
        self._clean_pdb_structure(lig_path)
        # now run ACPYPE on the ligand to produce the topology file
        self._parametrisation_pipeline(lig_path)

        # produces MOL.itp, need to separate the atomtypes directive out into ffMOL.itp for pmx
        # to generate the forcefield later
        self._separate_atomtypes(lig_path)

        # if we get through to here, return exit status 0

        q[node.get_node_id()] = 0
