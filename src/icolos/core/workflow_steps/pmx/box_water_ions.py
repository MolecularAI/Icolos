from typing import Dict, List
from icolos.core.containers.perturbation_map import Edge
from icolos.core.workflow_steps.pmx.base import StepPMXBase
from pydantic import BaseModel
import os
from icolos.utils.enums.program_parameters import (
    GromacsEnum,
    StepPMXEnum,
)
from icolos.core.workflow_steps.step import _LE
from icolos.utils.execute_external.pmx import PMXExecutor
from icolos.utils.general.parallelization import SubtaskContainer

_PSE = StepPMXEnum()
_GE = GromacsEnum()


class StepPMXBoxWaterIons(StepPMXBase, BaseModel):
    """
    Take the prepared structure files and prepare the system,
    runs editconf, solvate, genion and grompp for each system
    to be simulated
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
            run_func=self.boxWaterIons,
            step_id="pmx boxWaterIons",
            result_checker=self._check_result,
        )

    def boxWaterIons(self, jobs: List[str]):
        mdp_path = os.path.join(self.work_dir, "input/mdp")

        for edge in jobs:
            outLigPath = self._get_specific_path(
                workPath=self.work_dir, edge=edge, wp="ligand"
            )
            outProtPath = self._get_specific_path(
                workPath=self.work_dir, edge=edge, wp="complex"
            )

            # box ligand
            inStr = "{0}/init.pdb".format(outLigPath)
            outStr = "{0}/box.pdb".format(outLigPath)
            editconf_args = [
                "-f",
                inStr,
                "-o",
                outStr,
                "-bt",
                self.boxshape,
                "-d",
                self.boxd,
            ]
            self._gromacs_executor.execute(
                command=_GE.EDITCONF, arguments=editconf_args, check=True
            )
            # box protein
            inStr = "{0}/init.pdb".format(outProtPath)
            outStr = "{0}/box.pdb".format(outProtPath)
            editconf_args = [
                "-f",
                inStr,
                "-o",
                outStr,
                "-bt",
                self.boxshape,
                "-d",
                self.boxd,
            ]
            self._gromacs_executor.execute(
                command=_GE.EDITCONF, arguments=editconf_args, check=True
            )
            # water ligand
            inStr = "{0}/box.pdb".format(outLigPath)
            outStr = "{0}/water.pdb".format(outLigPath)
            top = "{0}/topol.top".format(outLigPath)
            solvate_args = ["-cp", inStr, "-cs", "spc216.gro", "-p", top, "-o", outStr]
            self._gromacs_executor.execute(command=_GE.SOLVATE, arguments=solvate_args)
            # water protein
            inStr = "{0}/box.pdb".format(outProtPath)
            outStr = "{0}/water.pdb".format(outProtPath)
            top = "{0}/topol.top".format(outProtPath)
            solvate_args = ["-cp", inStr, "-cs", "spc216.gro", "-p", top, "-o", outStr]
            self._gromacs_executor.execute(command=_GE.SOLVATE, arguments=solvate_args)

            # ions ligand
            inStr = "{0}/water.pdb".format(outLigPath)
            outStr = "{0}/ions.pdb".format(outLigPath)
            mdp = "{0}/em_l0.mdp".format(mdp_path)
            tpr = "{0}/tpr.tpr".format(outLigPath)
            top = "{0}/topol.top".format(outLigPath)
            mdout = "{0}/mdout.mdp".format(outLigPath)
            grompp_args = [
                "-f",
                mdp,
                "-c",
                inStr,
                "-r",
                inStr,
                "-p",
                top,
                "-o",
                tpr,
                "-maxwarn",
                4,
                "-po",
                mdout,
            ]

            self._gromacs_executor.execute(
                command=_GE.GROMPP, arguments=grompp_args, check=True
            )
            genion_args = [
                "-s",
                tpr,
                "-p",
                top,
                "-o",
                outStr,
                "-conc",
                self.conc,
                "-pname",
                self.pname,
                "-nname",
                self.nname,
                "-neutral",
            ]
            self._gromacs_executor.execute(
                command=_GE.GENION,
                arguments=genion_args,
                check=True,
                pipe_input="echo SOL",
            )
            # ions protein
            inStr = "{0}/water.pdb".format(outProtPath)
            outStr = "{0}/ions.pdb".format(outProtPath)
            mdp = "{0}/em_l0.mdp".format(mdp_path)
            tpr = "{0}/tpr.tpr".format(outProtPath)
            top = "{0}/topol.top".format(outProtPath)
            mdout = "{0}/mdout.mdp".format(outProtPath)
            grompp_args = [
                "-f",
                mdp,
                "-c",
                inStr,
                "-r",
                inStr,
                "-p",
                top,
                "-o",
                tpr,
                "-maxwarn",
                4,
                "-po",
                mdout,
            ]

            self._gromacs_executor.execute(
                command=_GE.GROMPP, arguments=grompp_args, check=True
            )
            genion_args = [
                "-s",
                tpr,
                "-p",
                top,
                "-o",
                outStr,
                "-conc",
                self.conc,
                "-pname",
                self.pname,
                "-nname",
                self.nname,
                "-neutral",
            ]
            self._gromacs_executor.execute(
                command=_GE.GENION,
                arguments=genion_args,
                check=True,
                pipe_input="echo SOL",
            )

            # clean backed files
            self._clean_backup_files(outLigPath)

            self._clean_backup_files(outProtPath)

    def _check_result(self, batch: List[List[str]]) -> List[List[bool]]:
        """
        Look in each hybridStrTop dir and check the output pdb files exist for the edges
        """
        output_files = [
            "ligand/tpr.tpr",
            "complex/tpr.tpr",
        ]
        results = []
        for subjob in batch:
            subjob_results = []
            for job in subjob:
                subjob_results.append(
                    all(
                        [
                            os.path.isfile(os.path.join(self.work_dir, job, f))
                            for f in output_files
                        ]
                    )
                )
            results.append(subjob_results)
        return results
