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
)

_PE = PMXEnum()
_PAE = PMXAtomMappingEnum()
_GE = GromacsEnum()
_SGE = StepGromacsEnum()


class StepPMXabfe(StepPMXBase, BaseModel):
    """Setup files for an ABFE calculation."""

    _gromacs_executor: GromacsExecutor = GromacsExecutor()

    def __init__(self, **data):
        super().__init__(**data)
        self._initialize_backend(PMXExecutor)
        self._check_backend_availability()
        self._gromacs_executor = GromacsExecutor(prefix_execution=_SGE.GROMACS_LOAD)

    def _separate_protein_ligand(self):
        # separate out protein and ligand lines from the written complex.pdb

        with open(os.path.join(self.work_dir, "complex.pdb"), "r") as f:
            lines = f.readlines()
        protein_lines = []
        ligand_lines = []
        # TODO: tighten up the logic for identifying the ligand here
        for line in lines:
            if "ATOM" in line:
                protein_lines.append(line)
            elif "HETATM" in line and "HOH" not in line:
                ligand_lines.append(line)

        with open(os.path.join(self.work_dir, "protein.pdb"), "w") as f:
            f.writelines(protein_lines)

        with open(os.path.join(self.work_dir, "MOL.pdb"), "w") as f:
            f.writelines(ligand_lines)

        os.remove(os.path.join(self.work_dir, "complex.pdb"))

    def execute(self):
        """
        Required inputs:
        + Protein.top, protein.gro
        + ligand.itp, ligand.gro


        Execution:
            - Separete protein and ligand from complex
            - run pdb2gmx on protein -> generate protein.top, ligand.grp
            - run acpype on ligand -> generate ligand.itp, ligand.gro
            - run pmx abfe to set up the system, done!
        """
        # use the same single dir setup as for the rest of the pmx pipeline

        assert self.work_dir is not None and os.path.isdir(self.work_dir)

        complex_file = self.data.generic.get_argument_by_extension(
            "pdb", rtn_file_object=True
        )
        complex_file.write(os.path.join(self.work_dir, "complex.pdb"), join=False)

        self._separate_protein_ligand()

        # parametrise the ligand, generate the itp files, top and gro files for the ligand
        self._parametrisation_pipeline(
            self.work_dir, include_gro=True, include_top=True
        )

        # parametrise protein
        self._parametrise_protein(protein="protein.pdb", path="", output="protein.gro")

        # run abfe

        args = {
            "-pt": "topol.top",
            "-lt": "MOL.itp",
            "-pc": "protein.gro",
            "-lc": "MOL_GMX.gro",
        }
        self._backend_executor.execute(
            command=_PE.ABFE,
            arguments=self.get_arguments(args),
            location=self.work_dir,
            check=True,
        )


help_string = """
pmx abfe -h
usage: pmx [-h] [-pt protop] [-lt ligtop] [-pc procrd] [-lc ligcrd] [--build]
           [--doublebox] [--longest_axis] [--keep_intra] [--lig_ids  ]
           [--pro_ids  ] [--restr_switch_on] [--seed int]

This scripts helps to setup an absolute binding free energy calculation. As a
minimal input, you need to provide a structure and topology file for both the
protein (or host) and ligand (or guest) molecule. The topology is setup so to
contain restraints as defined by Boresch et al. (2003) J Phys Chem B 107(35);
these include one distance, two angles, and three dihedrals between ligand and
protein. You can either provide explicitly the atoms to be included in the
restraints, or let the script choose them automatically.

optional arguments:
  -h, --help         show this help message and exit
  -pt protop         Input topology file for the protein. Default is
                     "protein.top".
  -lt ligtop         Input topology file for the ligand. It is expected that
                     all params needed for the ligand are explicitly defined
                     in this file. Default is "ligand.itp".
  -pc procrd         Input structure file in PDB or GRO format for the
                     protein. Default is "protein.gro".
  -lc ligcrd         Input structure file in PDB or GRO format for the ligand.
                     Default is "ligand.gro".
  --build            Whether to build the system (editconf, solvate, genion)
                     with a standard setup once the input files (top, gro) are
                     ready.
  --doublebox        Whether to use the double-system single-box setup. This
                     is useful for charged ligands. Default is False.
  --longest_axis     Whether to just place structures along the longest axis,
                     rather then minimising the volume. This option is
                     relevant only when using --doublebox. Default is False.
  --keep_intra       Whether to keep the LJ intramolecular interactions when
                     the ligand is decoupled. This option is relevant only
                     when using --doublebox. Default is False.
  --lig_ids          Three atom indices. If provided, these will be used for
                     the protein-ligand restraints. Otherwise they are chosen
                     automatically.
  --pro_ids          Three atom indices. If provided, these will be used for
                     the protein-ligand restraints. Otherwise they are chosen
                     automatically.
  --restr_switch_on  Whether to switch the restraints on or off, where "on"
                     means no restraints in stateA, and "off" means no
                     restraints in state B. Default is True (switch on).
  --seed int         Random seed to use when picking atoms for the restraints.
                     The automated restraints selection is stochastic, so if
                     you want to have a reproducible behaviour, provide a
                     random seed.
"""
