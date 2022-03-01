import shutil
from typing import List
from icolos.core.containers.compound import Compound
from icolos.core.workflow_steps.pmx.base import StepPMXBase
from pydantic import BaseModel
import os
from icolos.utils.enums.step_enums import StepGromacsEnum
from icolos.utils.execute_external.gromacs import GromacsExecutor
from icolos.utils.execute_external.pmx import PMXExecutor
from icolos.utils.enums.program_parameters import (
    GromacsEnum,
    PMXEnum,
    PMXAtomMappingEnum,
    StepPMXEnum,
)
from icolos.utils.general.parallelization import SubtaskContainer

_PE = PMXEnum()
_GE = GromacsEnum()
_SGE = StepGromacsEnum()
_PSE = StepPMXEnum()


class StepPMXabfe(StepPMXBase, BaseModel):
    """Setup files for an ABFE calculation."""

    _gromacs_executor: GromacsExecutor = GromacsExecutor()

    def __init__(self, **data):
        super().__init__(**data)
        self._initialize_backend(PMXExecutor)
        self._check_backend_availability()
        self._gromacs_executor = GromacsExecutor(prefix_execution=_SGE.GROMACS_LOAD)

    def execute(self):
        """
        This step manages the setup of a pmx ABFE run for a set of compounds and a protein target.
        Expects:
            - docked compounds provided as an sdf file
            - protein apo structure
            - directory containing mdp files
        """

        assert self.work_dir is not None and os.path.isdir(self.work_dir)
        replicas = (
            self.settings.additional["replicas"]
            if "replicas" in self.settings.additional.keys()
            else 3
        )
        # miror the dir structure for the input files used by the rbfe workflow
        os.makedirs(os.path.join(self.work_dir, "input"), exist_ok=True)
        for folder in ["ligands", "mdp", "protein"]:
            os.makedirs(os.path.join(self.work_dir, "input", folder), exist_ok=True)

        mdp_dir = self.data.generic.get_argument_by_extension(
            ext="mdp", rtn_file_object=True
        )
        # write mdp files to the input dir
        mdp_dir.write(os.path.join(self.work_dir, "input/mdp"))

        # create directory structure
        for comp in self.get_compounds():
            ident = comp.get_index_string()
            os.makedirs(os.path.join(self.work_dir, ident), exist_ok=True)

        # load in the provided apo structure
        protein = self.data.generic.get_argument_by_extension(
            "pdb", rtn_file_object=True
        )
        protein.write(
            os.path.join(self.work_dir, "input/protein/protein.pdb"), join=False
        )
        # parametrise protein
        self._parametrise_protein(
            protein="protein.pdb", path="input/protein", output="protein.gro"
        )
        # self._clean_protein()

        self.execution.parallelization.max_length_sublists = 1
        self._subtask_container = SubtaskContainer(
            max_tries=self.execution.failure_policy.n_tries
        )
        self._subtask_container.load_data(self.get_compounds())
        self._execute_pmx_step_parallel(
            run_func=self._parametrise_nodes,
            step_id="parametrize ligands",
            result_checker=self._check_params,
        )

        # now
        self._subtask_container.load_data(self.get_compounds())
        self._execute_pmx_step_parallel(
            run_func=self._setup_abfe,
            step_id="pxm abfe",
            result_checker=self._find_nan_vals,
        )

        # now we make the rest of the dir structure
        for comp in self.get_compounds():
            ident = comp.get_index_string()
            shutil.copyfile(
                os.path.join(self.work_dir, "input/protein/posre.itp"),
                os.path.join(self.work_dir, ident, "complex/posre.itp"),
            )
            for wp in self.therm_cycle_branches:
                wppath = os.path.join(self.work_dir, ident, wp)
                # copy the posre.itp file from input/protein to each

                # stateA/stateB - coupled + decoupled states
                for state in self.states:
                    statepath = os.path.join(wppath, state)
                    os.makedirs(statepath, exist_ok=True)

                    # run1/run2/run3
                    for r in range(1, replicas + 1):
                        runpath = os.path.join(statepath, f"run{r}")
                        os.makedirs(runpath, exist_ok=True)

                        # em/eq_posre/eq/transitions
                        # TODO: this differs from the equil used for RBFE - can we get away without the extra equilibration?
                        for sim in self.settings.additional[_PSE.SIM_TYPES]:
                            simpath = os.path.join(runpath, sim)
                            os.makedirs(simpath, exist_ok=True)

    def _setup_abfe(self, jobs):
        """
        Executes pmx abfe, moves the resulting built files to the right dir
        """
        if isinstance(jobs, list):
            comp = jobs[0]

        args = {
            "-pt": os.path.join(self.work_dir, "input/protein/topol.top"),
            "-lt": os.path.join(
                self.work_dir,
                "input/ligands",
                comp.get_index_string(),
                "MOL.acpype/MOL_GMX.itp",
            ),
            "-pc": os.path.join(self.work_dir, "input/protein/protein.gro"),
            "-lc": os.path.join(
                self.work_dir,
                "input/ligands",
                comp.get_index_string(),
                "MOL.acpype/MOL_GMX.gro",
            ),
            "--build": "",
        }
        self._backend_executor.execute(
            command=_PE.ABFE,
            arguments=self.get_arguments(args),
            location=os.path.join(self.work_dir, comp.get_index_string()),
            check=True,
        )
        # note that this is stochastic, and sometimes it generates bad restraaints/nan values
        # we will simply resubmit n times

    def _find_nan_vals(self, next_batch: List[str]) -> List[List[bool]]:
        """
        Looks throuh the dirs specified in jobs, reads restraints.info
        """
        # sublist length set to 1 for this step
        batch_results = []
        for subtask in next_batch:
            subtask_results = []
            for comp in subtask:
                with open(
                    os.path.join(
                        self.work_dir, comp.get_index_string(), "restraints.info"
                    ),
                    "r",
                ) as f:
                    lines = f.readlines()

                subtask_results.append(any(["nan" in l for l in lines]))
            batch_results.append(subtask_results)
        return batch_results

    def _check_params(self, batch: List[List[Compound]]) -> List[List[bool]]:
        """
        check ligand parameters have been generated properly
        """
        output_files = ["MOL.itp", "MOL.mol2"]
        results = []
        for subjob in batch:
            subjob_results = []
            for job in subjob:
                subjob_results.append(
                    all(
                        [
                            os.path.isfile(
                                os.path.join(self.work_dir, job.get_index_string(), f)
                            )
                            for f in output_files
                        ]
                    )
                )
            results.append(subjob_results)
        return results
