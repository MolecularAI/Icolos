from rdkit import Chem


def to_smiles(mol, isomericSmiles=False):
    """
    Converts a Mol object into a canonical SMILES string.
    :param mol: Mol object.
    :return: A SMILES string.
    """
    return Chem.MolToSmiles(mol, isomericSmiles=isomericSmiles)
