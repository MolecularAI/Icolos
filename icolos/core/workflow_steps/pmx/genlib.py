from icolos.core.workflow_steps.pmx.base import StepPMXBase
from pydantic import BaseModel


class StepPMXgenlib(StepPMXBase, BaseModel):
    """Generate pmx ff library."""

    def __init__(self, **data):
        super().__init__(**data)

    def execute(self):
        pass


help_string = """
pmx genlib -h
usage: pmx [-h] [-f1 ipdb1] [-f2 ipdb2] [-o1 opdb1] [-o2 opdb2]
           [--ffpath ffpath] [--fatp fatp] [--fnb fnb] [--moltype moltype]
           [--noalign] [--cbeta] [--noH2Heavy] [--log log]

The script creates hybrid structure and topology database entries (mtp and rtp)
in order to generate a pmx alchemical force field library.

The easiest way to generate the library is to call this script from within
the folder of the force field you are interested in.

If two pdb files (aligned on the backbone) are provided, the hybrid pdb, mtp,
and rtp files are written to file. If no pdb input file is provided,
the script uses pregenerated residues in order to build hybrid pdb, mtp, and
rtp files for all possible residue pairs, thus preparing the whole pmx ff
library.

In addition, atomtype (-fatp) and non-bonded parameter (-fnm) files for the
introduced dummy atoms are generated. By default, these point towards the
files already present in the forcefield. In this way, the additional parameters
for the dummies are appended to the existing ff file, rather than being
written to new files.

optional arguments:
  -h, --help         show this help message and exit
  -f1 ipdb1          First input PDB file. Default is none provided.
  -f2 ipdb2          Second input PDB file. Default is none provided.
  -o1 opdb1          First output PDB file. Default is none provided.
  -o2 opdb2          Second output PDB file. Default is none provided.
  --ffpath ffpath    Path to mutation forcefield. Default is current folder.
  --fatp fatp        Atom types (atp) file. If the file is 
                     present, data is appended to it, otherwise a new 
                     file is created. Default is "atomtypes.atp".
  --fnb fnb          Non-bonded (nb) types file. If the file is 
                     present, data is appended to it, otherwise a new 
                     file is created. Default is "ffnonbonded.itp".
  --moltype moltype  The type of molecule for which the library is 
                     being built. Available options are "protein", "dna", 
                     or "rna". Default is "protein".
  --noalign          Whether to align the sidechains of the two 
                     input PDB files provided. Default it True; this flag 
                     sets it to False.
  --cbeta            Whether to morph sidechain between the two 
                     residues or to use dummy atoms to (de)couple the 
                     whole sidechain. By default, sidechain atoms are 
                     morphed so to minimise the size of the perturbation. 
                     With this flag set, whole sidechains are (de)coupled 
                     instead; i.e. all atoms after C-beta are not mapped 
                     between the two residues.
  --noH2Heavy        Whether to allow hydrogen to/from heavy atoms 
                     morphing. Default is True, this flag sets it to False.
  --log log          Logging level. Either "info" or "debug". Default is "info".
"""
