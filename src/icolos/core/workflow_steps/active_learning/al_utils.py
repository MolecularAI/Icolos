import os
import shutil
import tempfile
from typing import List
import numpy as np

# RDkit
from rdkit import Chem
from rdkit.Chem.rdmolops import GetAdjacencyMatrix
from ase import io
import ase
from dscribe.descriptors import SOAP

# Pytorch and Pytorch Geometric
import torch
from torch_geometric.data import Data
from torch.utils.data import DataLoader

# functions adapted from https://www.blopig.com/blog/2022/02/how-to-turn-a-smiles-string-into-a-molecular-graph-for-pytorch-geometric/


def one_hot_encoding(x, permitted_list):
    """
    Maps input elements x which are not in the permitted list to the last element
    of the permitted list.
    """
    if x not in permitted_list:
        x = permitted_list[-1]
    binary_encoding = [
        int(boolean_value)
        for boolean_value in list(map(lambda s: x == s, permitted_list))
    ]
    return binary_encoding


def get_atom_features(atom, use_chirality=True, hydrogens_implicit=True):
    """
    Takes an RDKit atom object as input and gives a 1d-numpy array of atom features as output.
    """
    # define list of permitted atoms

    permitted_list_of_atoms = [
        "C",
        "N",
        "O",
        "S",
        "F",
        "Si",
        "P",
        "Cl",
        "Br",
        "Mg",
        "Na",
        "Ca",
        "Fe",
        "As",
        "Al",
        "I",
        "B",
        "V",
        "K",
        "Tl",
        "Yb",
        "Sb",
        "Sn",
        "Ag",
        "Pd",
        "Co",
        "Se",
        "Ti",
        "Zn",
        "Li",
        "Ge",
        "Cu",
        "Au",
        "Ni",
        "Cd",
        "In",
        "Mn",
        "Zr",
        "Cr",
        "Pt",
        "Hg",
        "Pb",
        "Unknown",
    ]

    if hydrogens_implicit == False:
        permitted_list_of_atoms = ["H"] + permitted_list_of_atoms

    # compute atom features

    atom_type_enc = one_hot_encoding(str(atom.GetSymbol()), permitted_list_of_atoms)

    n_heavy_neighbors_enc = one_hot_encoding(
        int(atom.GetDegree()), [0, 1, 2, 3, 4, "MoreThanFour"]
    )

    formal_charge_enc = one_hot_encoding(
        int(atom.GetFormalCharge()), [-3, -2, -1, 0, 1, 2, 3, "Extreme"]
    )

    hybridisation_type_enc = one_hot_encoding(
        str(atom.GetHybridization()),
        ["S", "SP", "SP2", "SP3", "SP3D", "SP3D2", "OTHER"],
    )

    is_in_a_ring_enc = [int(atom.IsInRing())]

    is_aromatic_enc = [int(atom.GetIsAromatic())]

    atomic_mass_scaled = [float((atom.GetMass() - 10.812) / 116.092)]

    vdw_radius_scaled = [
        float((Chem.GetPeriodicTable().GetRvdw(atom.GetAtomicNum()) - 1.5) / 0.6)
    ]

    covalent_radius_scaled = [
        float((Chem.GetPeriodicTable().GetRcovalent(atom.GetAtomicNum()) - 0.64) / 0.76)
    ]
    atom_feature_vector = (
        atom_type_enc
        + n_heavy_neighbors_enc
        + formal_charge_enc
        + hybridisation_type_enc
        + is_in_a_ring_enc
        + is_aromatic_enc
        + atomic_mass_scaled
        + vdw_radius_scaled
        + covalent_radius_scaled
    )

    if use_chirality == True:
        chirality_type_enc = one_hot_encoding(
            str(atom.GetChiralTag()),
            [
                "CHI_UNSPECIFIED",
                "CHI_TETRAHEDRAL_CW",
                "CHI_TETRAHEDRAL_CCW",
                "CHI_OTHER",
            ],
        )
        atom_feature_vector += chirality_type_enc

    if hydrogens_implicit == True:
        n_hydrogens_enc = one_hot_encoding(
            int(atom.GetTotalNumHs()), [0, 1, 2, 3, 4, "MoreThanFour"]
        )
        atom_feature_vector += n_hydrogens_enc
    return np.array(atom_feature_vector)


def get_bond_features(bond, use_stereochemistry=True):
    """
    Takes an RDKit bond object as input and gives a 1d-numpy array of bond features as output.
    """
    permitted_list_of_bond_types = [
        Chem.rdchem.BondType.SINGLE,
        Chem.rdchem.BondType.DOUBLE,
        Chem.rdchem.BondType.TRIPLE,
        Chem.rdchem.BondType.AROMATIC,
    ]
    bond_type_enc = one_hot_encoding(bond.GetBondType(), permitted_list_of_bond_types)

    bond_is_conj_enc = [int(bond.GetIsConjugated())]

    bond_is_in_ring_enc = [int(bond.IsInRing())]

    bond_feature_vector = bond_type_enc + bond_is_conj_enc + bond_is_in_ring_enc

    if use_stereochemistry == True:
        stereo_type_enc = one_hot_encoding(
            str(bond.GetStereo()), ["STEREOZ", "STEREOE", "STEREOANY", "STEREONONE"]
        )
        bond_feature_vector += stereo_type_enc
    return np.array(bond_feature_vector)


def create_graph(mol: Chem.Mol, y):
    """
    Inputs:

    x_smiles = [smiles_1, smiles_2, ....] ... a list of SMILES strings
    y = [y_1, y_2, ...] ... a list of numerial labels for the SMILES strings (such as associated pKi values)

    Outputs:

    data_list = [G_1, G_2, ...] ... a list of torch_geometric.data.Data objects which represent labeled molecular graphs that can readily be used for machine learning

    """

    data_list = []

    # convert SMILES to RDKit mol object

    # get feature dimensions
    n_nodes = mol.GetNumAtoms()
    n_edges = 2 * mol.GetNumBonds()
    unrelated_smiles = "O=O"
    unrelated_mol = Chem.MolFromSmiles(unrelated_smiles)
    n_node_features = len(get_atom_features(unrelated_mol.GetAtomWithIdx(0)))
    n_edge_features = len(get_bond_features(unrelated_mol.GetBondBetweenAtoms(0, 1)))
    # construct node feature matrix X of shape (n_nodes, n_node_features)
    X = np.zeros((n_nodes, n_node_features))
    for atom in mol.GetAtoms():
        X[atom.GetIdx(), :] = get_atom_features(atom)

    X = torch.tensor(X, dtype=torch.float)

    # construct edge index array E of shape (2, n_edges)
    (rows, cols) = np.nonzero(GetAdjacencyMatrix(mol))
    torch_rows = torch.from_numpy(rows.astype(np.int64)).to(torch.long)
    torch_cols = torch.from_numpy(cols.astype(np.int64)).to(torch.long)
    E = torch.stack([torch_rows, torch_cols], dim=0)

    # construct edge feature array EF of shape (n_edges, n_edge_features)
    EF = np.zeros((n_edges, n_edge_features))

    for (k, (i, j)) in enumerate(zip(rows, cols)):

        EF[k] = get_bond_features(mol.GetBondBetweenAtoms(int(i), int(j)))

    EF = torch.tensor(EF, dtype=torch.float)

    # construct label tensor
    y_tensor = torch.tensor(np.array([y]), dtype=torch.float)

    # construct Pytorch Geometric data object and append to data list
    data_list.append(Data(x=X, edge_index=E, edge_attr=EF, y=y_tensor))
    return DataLoader(data_list, batch_size=32)


def greedy_acquisition(
    estimator,
    X: np.ndarray,
    previous_idx,
    n_instances: int,
) -> np.ndarray:
    """
    Implement greedy acquisition strategy, return the n_samples best scores

    """
    try:
        predictions = estimator.predict(X)
    except:
        print("estimator not fitted, producing random ")
        predictions = np.random.uniform(12, 0, X.shape[0])

    # zero those predictions we've seen before
    for idx in previous_idx:
        predictions[idx] = 0
    # smaller before n_instances, largest after
    sorted_preds = np.argpartition(predictions, -n_instances, axis=0)[-n_instances:]
    return sorted_preds


def convert_mol_to_ase_atoms(mol: Chem.Mol) -> ase.Atoms:
    tmp_dir = tempfile.mkdtemp()
    Chem.rdmolfiles.MolToXYZFile(mol, os.path.join(tmp_dir, "mol.xyz"))
    atoms = io.read(os.path.join(tmp_dir, "mol.xyz"))
    shutil.rmtree(tmp_dir)
    return atoms
