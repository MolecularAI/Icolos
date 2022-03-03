import os
import gpytorch
from dscribe.kernels import REMatchKernel
import torch
from rdkit import Chem
from dscribe.descriptors import SOAP
from sklearn.preprocessing import normalize

from icolos.core.workflow_steps.active_learning.al_utils import (
    convert_mol_to_ase_atoms,
)

# define a custom kernel to compute similarity between two compounds based on their soap matrices


class SOAP_Kernel(gpytorch.kernels.Kernel):
    def __init__(self, tmp_dir: str, device: str) -> None:
        super().__init__(active_dims=[0])
        self.tmp_dir = tmp_dir
        self.device = device
        is_stationary = True

    def forward(self, x1: torch.Tensor, x2: torch.Tensor = None, **kwargs):

        atoms = []
        for i in x1.flatten():
            mol = Chem.SDMolSupplier(
                os.path.join(self.tmp_dir, f"{i}.sdf"), removeHs=True
            )[0]
            atoms.append(convert_mol_to_ase_atoms(mol))
        soap = SOAP(
            species=["C", "H", "O", "N", "F", "Cl", "I", "Br", "S", "P"],
            rcut=5,
            nmax=8,
            lmax=6,
            crossover=True,
            sparse=True,
            dtype="float32",
        )

        soap_descriptors = []
        for sys in atoms:
            desc = soap.create(sys).to_scipy_sparse()
            soap_descriptors.append(normalize(desc, copy=False, axis=1))
        # TODO: this kernel computation needs to be parallelized
        soap_k = REMatchKernel(
            # TODO: integration is much slower than the gaussian radial functions, switch if performance becomes an issue
            metric="polynomial",
            degree=3,
            gamma=1,
            coef0=0,
            alpha=0.5,
            threshold=1e-6,
            normalize_kernel=True,
        )
        # TODO: augment this with something to conver smoothness, Matern kernel or something
        return torch.from_numpy(soap_k.create(soap_descriptors)).to(self.device)


class SOAP_GP(gpytorch.models.ExactGP):
    def __init__(self, likelihood, tmp_dir: str, device: str, noise_init=None):

        super().__init__(train_inputs=None, train_targets=None, likelihood=likelihood)
        self.mean_module = gpytorch.means.ConstantMean()
        self.covar_module = SOAP_Kernel(tmp_dir=tmp_dir, device=device)
        self.device = device

    def forward(self, x, *args, **kwargs):
        # x is a matrix of identifiers for the compounds being passed, let the kernel compute these on the fly, should of of shape (n_compounds, 1), where the last dim is the identifier or something
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)
