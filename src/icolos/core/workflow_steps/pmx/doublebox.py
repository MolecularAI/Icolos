from icolos.core.workflow_steps.pmx.base import StepPMXBase
from pydantic import BaseModel


class StepPMXdoublebox(StepPMXBase, BaseModel):
    """Place two input structures into a single box."""

    def __init__(self, **data):
        super().__init__(**data)

    def execute(self):
        pass


help_string = """
pmx doublebox -h
usage: pmx [-h] -f1  -f2  [-o] [-r] [-d] [--longest_axis]

Places two structures into a single box. The box is a rectangular cuboid in
which the two structures are placed in such a way as to minimise the box
volume. You can use this script to help in the setup of a calculation using
the single-box double-system approach.

optional arguments:
  -h, --help      show this help message and exit
  -f1             First structure in PDB or GRO format.
  -f2             Second structure in PDB or GRO format.
  -o              Name of output file. Default is "doublebox.pdb".
  -r              Distance between the two structures (nm). Default is 2.5 nm.
  -d              Distance to the box wall (nm). Default is 1.5 nm.
  --longest_axis  Whether to just place structures along the longest axis,
                  rather then minimising the volume. Default is False.
"""
