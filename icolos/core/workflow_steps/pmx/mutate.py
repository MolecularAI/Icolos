from icolos.core.workflow_steps.pmx.base import StepPMXBase
from pydantic import BaseModel


class StepPMXmutate(StepPMXBase, BaseModel):
    """Mutate protein or DNA/RNA."""

    def __init__(self, **data):
        super().__init__(**data)

    def execute(self):
        pass


help_string = """
pmx mutate -h
usage: pmx [-h] [-f infile] [-fB infileB] [-o outfile] [-ff ff]                                    
           [--script script] [--keep_resid | --ref ] [--resinfo]                                   
                                                                                                   
This script applies mutations of residues in a structure file for subsequent                                          
free energy calculations. It supports mutations to protein, DNA, and RNA                           
molecules.                                                                                         
                                                                          
The mutation information and dummy placements are taken from the hybrid residue
database "mutres.mtp". The best way to use this script is to take a pdb/gro file                                                                    
that has been written with pdb2gmx with all hydrogen atoms present.
                                                                          
By default, all residues are renumbered starting from 1, so to have unique
residue IDs. If you want to keep the original residue IDs, you can use the flag                                                                     
--keep_resid. In this case, you will also need to provide chain information                                                                         
in order to be able to mutate the desired residue. Alternatively, if you would                                                                      
like to use the original residue IDs but these have been changed, e.g. by gromacs,
you can provide a reference PDB file (with chain information too) using the --ref
flag. The input structure will be mutated according to the IDs chosen for the
reference structure after having mapped the two residue indices.

The program can either be executed interactively or via script. The script file
simply has to consist of "residue_id target_residue_name" pairs (just with some
space between the id and the name), or "chain_id residue_id target_residue_name"
if you are keeping the original residue IDs or providing a reference structure.

The script uses an extended one-letter code for amino acids to account for
different protonation states. Use the --resinfo flag to print the dictionary.

optional arguments:
  -h, --help       show this help message and exit
  -f infile        Input structure file in PDB or GRO format. Default is "protein.pdb"
  -fB infileB      Input structure file of the B state in PDB or GRO format (optional).
  -o outfile       Output structure file in PDB or GRO format. Default is "mutant.pdb"
  -ff ff           Force field to use. If none is provided, 
                   a list of available ff will be shown.
  --script script  Text file with list of mutations (optional).
  --keep_resid     Whether to renumber all residues or to keep the
                   original residue IDs. By default, all residues are
                   renumbered so to have unique IDs. With this flags set,
                   the original IDs are kept. Because the IDs might not
                   be unique anymore, you will also be asked to choose
                   the chain ID where the residue you want to mutate is.
  --ref            Provide a reference PDB structure from which to map
                   the chain and residue IDs onto the file to be mutated (-f).
                   This can be useful when wanting to mutate a file that
                   has had its residues renumbered or the chain information
                   removed (e.g. after gmx grompp). As in the --keep_resid
                   option, if --ref is chosen, you will need to provide chain
                   information either interactively or via the --script flag.
  --resinfo        Show the list of 3-letter -> 1-letter residues
"""
