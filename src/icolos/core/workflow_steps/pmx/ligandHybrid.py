from typing import List
from icolos.core.containers.perturbation_map import Edge
from icolos.core.workflow_steps.pmx.base import StepPMXBase
from pydantic import BaseModel
import os
from icolos.utils.enums.program_parameters import PMXEnum, PMXAtomMappingEnum
from icolos.core.workflow_steps.step import _LE
import numpy as np

from icolos.utils.execute_external.pmx import PMXExecutor
from icolos.utils.general.parallelization import SubtaskContainer

_PE = PMXEnum()
_PAE = PMXAtomMappingEnum()


class StepPMXligandHybrid(StepPMXBase, BaseModel):
    """Ligand alchemy: hybrid structure/topology."""

    def __init__(self, **data):
        super().__init__(**data)
        self._initialize_backend(PMXExecutor)

    def _execute_command(self, jobs: List[str]):
        edge = jobs[0]
        parts = edge.split("_")
        lig1 = parts[0]
        lig2 = parts[1]

        arguments = {
            "-i1": os.path.join(
                self.work_dir,
                _PAE.LIGAND_DIR,
                lig1,
                "MOL.pdb",
            ),
            "-i2": os.path.join(
                self.work_dir,
                _PAE.LIGAND_DIR,
                lig2,
                "MOL.pdb",
            ),
            "-itp1": os.path.join(self.work_dir, _PAE.LIGAND_DIR, lig1, "MOL.itp"),
            "-itp2": os.path.join(self.work_dir, _PAE.LIGAND_DIR, lig2, "MOL.itp"),
        }
        # write output files the hybridStrTop directory for each edge
        output_dir = os.path.join(self.work_dir, edge, _PE.HYBRID_STR_TOP)
        arguments = self._prepare_arguments(args=arguments, output_dir=output_dir)
        self._backend_executor.execute(
            command=_PE.LIGANDHYBRID,
            arguments=arguments,
            check=True,
            location=self.work_dir,
        )

    def _prepare_arguments(self, args, output_dir):
        """
        Prepare the final set of arguments as a list, config overrides defaults
        """
        prepared_args = []
        default_args = {
            "-pairs": f"{output_dir}/pairs1.dat",
            "-oA": f"{output_dir}/mergedA.pdb",
            "-oB": f"{output_dir}/mergedB.pdb",
            "-oitp": f"{output_dir}/merged.itp",
            "-offitp": f"{output_dir}/ffmerged.itp",
            "-log": f"{output_dir}/mapping.log",
        }
        for key, value in args.items():
            default_args[key] = value

        for key, value in self.settings.arguments.parameters.items():
            default_args[key] = value

        for key, value in default_args.items():
            prepared_args.append(key),
            prepared_args.append(value)

        for flag in self.settings.arguments.flags:
            prepared_args.append(flag)
        return prepared_args

    def execute(self):
        assert self.work_dir is not None and os.path.isdir(self.work_dir)

        edges = [e.get_edge_id() for e in self.get_edges()]

        self.execution.parallelization.max_length_sublists = 1
        self._subtask_container = SubtaskContainer(
            max_tries=self.execution.failure_policy.n_tries
        )
        self._subtask_container.load_data(edges)
        self._execute_pmx_step_parallel(
            run_func=self._execute_command,
            step_id="ligandHybrid",
            result_checker=self._check_result,
        )

    def _check_result(self, batch: List[List[str]]) -> List[List[bool]]:
        """
        Look in each hybridStrTop dir and check the output pdb files exist for the edges
        """
        output_files = ["ffmerged.itp", "mergedB.pdb", "mergedA.pdb"]
        results = []
        for subjob in batch:
            subjob_results = []
            for job in subjob:
                subjob_results.append(
                    all(
                        [
                            os.path.isfile(
                                os.path.join(self.work_dir, job, "hybridStrTop", f)
                            )
                            for f in output_files
                        ]
                    )
                )
            results.append(subjob_results)
        return results


help_string = """
pmx ligandHybrid -h
usage: pmx [-h] [-i1 lig1.pdb] [-i2 lig2.pdb] [-itp1 lig1.itp]
           [-itp2 lig2.itp] [-pairs pairs.dat] [-n1 scaffold1.ndx]
           [-n2 scaffold2.ndx] [-oA mergedA.pdb] [-oB mergedB.pdb]
           [-oitp merged.itp] [-offitp ffmerged.itp] [-log hybrid.log]
           [--d 0.05] [--fit] [--split] [--scDUMm 1.0] [--scDUMa 1.0]
           [--scDUMd 1.0] [--deAng]

Provided two structures and topologies, build hybrid structure/topology.

optional arguments:
  -h, --help            show this help message and exit
  -i1 lig1.pdb          Input ligand structure 1. Default is "lig1.pdb"
  -i2 lig2.pdb          Input ligand structure 2. Default is "lig2.pdb"
  -itp1 lig1.itp        Input ligand topology 1. Default is "lig1.itp"
  -itp2 lig2.itp        Input ligand topology 2. Default is "lig2.itp"
  -pairs pairs.dat      Optional input: atom pair mapping. 
  -n1 scaffold1.ndx     Optional input: index of atoms to consider for mol1 
  -n2 scaffold2.ndx     Optional input: index of atoms to consider for mol2 
  -oA mergedA.pdb       Output: hybrid structure based on the ligand 1. Default is "mergedA.pdb"
  -oB mergedB.pdb       Output: hybrid structure based on the ligand 2. Default is "mergedB.pdb"
  -oitp merged.itp      Output: hybrid topology. Default is "merged.itp"
  -offitp ffmerged.itp  Output: atomtypes for hybrid topology. Default is "ffmerged.itp"
  -log hybrid.log       Output: log file. Default is "hybrid.log"
  --d 0.05              Optional: if -pairs not provided, distance (nm) between atoms to consider them morphable for alignment approach (default 0.05 nm).
  --fit                 Fit mol2 onto mol1, only works if pairs.dat is provided
  --split               split the topology into separate transitions
  --scDUMm 1.0          scale dummy masses using the counterpart atoms
  --scDUMa 1.0          scale bonded dummy angle parameters
  --scDUMd 1.0          scale bonded dummy dihedral parameters
  --deAng               decouple angles composed of 1 dummy and 2 non-dummies
"""
