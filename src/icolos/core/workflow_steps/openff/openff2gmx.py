from icolos.core.workflow_steps.step import StepBase
from pydantic import BaseModel
from openmm.app import PDBFile

import parmed
from openff.toolkit.topology import Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.utils import get_data_file_path

from icolos.utils.enums.step_enums import StepOpenFFEnum
import os

# Note that this implementation leverages parmed for now, although will likely move over to the OpenFF Interchange tooling once stable
_SOFE = StepOpenFFEnum()

# TOOD: Eventually this step should be able to do a pdb2gmx job
class StepOFF2gmx(StepBase, BaseModel):
    def __init__(self, **data):
        super().__init__(**data)

    def parametrise_mols(self, tmp_dir: str):
        """
        Generate parameters for each mol
        """
        # TODO: do we want to throw everything together or split the params into separate files?
        mols = [
            Molecule.from_smiles(smi)
            for smi in self.get_additional_setting(_SOFE.UNIQUE_MOLS)
        ]
        pdb_file = PDBFile(
            os.path.join(tmp_dir, self.data.generic.get_argument_by_extension("pdb"))
        )
        omm_topology = pdb_file.topology

        off_topology = Topology.from_openmm(omm_topology, unique_molecules=mols)

        forcefield = ForceField(self.get_additional_setting(_SOFE.FORCEFIELD))

        omm_system = forcefield.create_openmm_system(off_topology)

        if self.get_additional_setting(_SOFE.METHOD, _SOFE.PARMED) == _SOFE.PARMED:
            parmed_struct = parmed.openmm.load_topology(
                omm_topology, omm_system, pdb_file.positions
            )
            parmed_struct.save(os.path.join(tmp_dir, "MOL.top"), overwrite=True)
            parmed_struct.save(os.path.join(tmp_dir, "MOL.gro"), overwrite=True)

            # TODO: validate energy differences

        elif self.get_additional_setting(_SOFE.METHOD) == _SOFE.INTERCHANGE:
            raise NotImplementedError
        else:
            raise NotImplementedError

    def execute(self):
        """
        Builds a system and parametrise using OpenFF SAGE params, then convert to a GROMACS top/gro format for downstream simulation
        """
        tmp_dir = self._make_tmpdir()
        self.data.generic.write_out_all_files(tmp_dir)

        self.parametrise_mols(tmp_dir)

        self._parse_output(tmp_dir)
        self._remove_temporary(tmp_dir)


# If we want to build OpenFF params instead of gaff, we would need to make a call to a different parametrisation pipeline, then load the protein into a ParmEd System, combine with the OpenFF-built system and combine into a gro/top file.
