import random

import rdkit.Chem as rkc
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import SaltRemover
from rdkit.Chem import rdmolops


def _initialiseNeutralisationReactions():
    patts = (
        # Imidazoles
        ("[n+;H]", "n"),
        # Amines
        ("[N+;!H0]", "N"),
        # Carboxylic acids and alcohols
        ("[$([O-]);!$([O-][#7])]", "O"),
        # Thiols
        ("[S-;X1]", "S"),
        # Sulfonamides
        ("[$([N-;X2]S(=O)=O)]", "N"),
        # Enamines
        ("[$([N-;X2][C,N]=C)]", "N"),
        # Tetrazoles
        ("[n-]", "[nH]"),
        # Sulfoxides
        ("[$([S-]=O)]", "S"),
        # Amides
        ("[$([N-]C=O)]", "N"),
    )
    return [(Chem.MolFromSmarts(x), Chem.MolFromSmiles(y, False)) for x, y in patts]


def _neutralise_charges(mol, reactions=None):
    if reactions is None:
        reactions = _initialiseNeutralisationReactions()
    replaced = False
    for i, (reactant, product) in enumerate(reactions):
        while mol.HasSubstructMatch(reactant):
            replaced = True
            rms = AllChem.ReplaceSubstructs(mol, reactant, product)
            mol = rms[0]
    if replaced:
        return mol, True
    else:
        return mol, False


def _get_largest_fragment(mol):
    frags = rdmolops.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    maxmol = None
    for mol in frags:
        if mol is None:
            continue
        if maxmol is None:
            maxmol = mol
        if maxmol.GetNumHeavyAtoms() < mol.GetNumHeavyAtoms():
            maxmol = mol
    return maxmol


_saltremover = SaltRemover.SaltRemover()


def _valid_size(
    mol, min_heavy_atoms, max_heavy_atoms, element_list, remove_long_side_chains
):
    """Filters molecules on number of heavy atoms and atom types"""
    mol = _rare_filters(mol)
    if mol:
        correct_size = min_heavy_atoms < mol.GetNumHeavyAtoms() < max_heavy_atoms
        if not correct_size:
            return

        valid_elements = all(
            [atom.GetAtomicNum() in element_list for atom in mol.GetAtoms()]
        )
        if not valid_elements:
            return

        has_long_sidechains = False
        if remove_long_side_chains:
            # remove aliphatic side chains with at least 5 carbons not in a ring
            sma = "[CR0]-[CR0]-[CR0]-[CR0]-[CR0]"
            has_long_sidechains = mol.HasSubstructMatch(Chem.MolFromSmarts(sma))

        return correct_size and valid_elements and not has_long_sidechains


def _rare_filters(mol):
    if mol:
        ciano_filter = "[C-]#[N+]"
        oh_filter = "[OH+]"
        sulfur_filter = "[SH]"
        if (
            not mol.HasSubstructMatch(Chem.MolFromSmarts(ciano_filter))
            and not mol.HasSubstructMatch(Chem.MolFromSmarts(oh_filter))
            and not mol.HasSubstructMatch(Chem.MolFromSmarts(sulfur_filter))
        ):
            return mol


def convert_to_rdkit_smiles(smiles):
    return Chem.MolToSmiles(
        Chem.MolFromSmiles(smiles, sanitize=False), isomericSmiles=True
    )


def randomize_smiles(smiles, random_type="restricted"):
    """
    Returns a random SMILES given a SMILES of a molecule.
    :param random_type: The type (unrestricted, restricted) of randomization performed.
    :return : A random SMILES string of the same molecule or None if the molecule is invalid.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None

    if random_type == "unrestricted":
        return rkc.MolToSmiles(mol, canonical=False, doRandom=True, isomericSmiles=True)
    if random_type == "restricted":
        new_atom_order = list(range(mol.GetNumHeavyAtoms()))
        random.shuffle(new_atom_order)
        random_mol = rkc.RenumberAtoms(mol, newOrder=new_atom_order)
        return rkc.MolToSmiles(random_mol, canonical=False, isomericSmiles=True)
    raise ValueError("Type '{}' is not valid".format(random_type))


def to_mol(smi):
    """
    Creates a Mol object from a SMILES string.
    :param smi: SMILES string.
    :return: A Mol object or None if it's not valid.
    """
    if smi:
        return rkc.MolFromSmiles(smi)


def to_smiles(mol):
    """
    Converts a Mol object into a canonical SMILES string.
    :param mol: Mol object.
    :return: A SMILES string.
    """
    return rkc.MolToSmiles(mol, isomericSmiles=True)
