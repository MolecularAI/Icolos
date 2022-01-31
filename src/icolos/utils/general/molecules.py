from rdkit import Chem
from rdkit.Chem import rdmolops

from icolos.utils.enums.compound_enums import CompoundTagsEnum


def get_charge_for_molecule(molecule: Chem.Mol, add_as_tag=False) -> int:
    _MTE = CompoundTagsEnum()
    charge = rdmolops.GetFormalCharge(molecule)
    if add_as_tag:
        molecule.SetProp(_MTE.FORMAL_CHARGE_TAG, str(charge))
    return charge


def write_molecule_to_sdf(path: str, molecule: Chem.Mol):
    if molecule is None or not isinstance(molecule, Chem.Mol):
        raise ValueError("Function requires input attribute to be an RDkit molecule.")
    writer = Chem.SDWriter(path)
    writer.write(molecule)
    writer.close()
