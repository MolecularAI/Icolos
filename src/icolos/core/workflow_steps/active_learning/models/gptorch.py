import gpytorch
from dscribe.kernels import REMatchKernel
import torch
import numpy as np
from sklearn.preprocessing import normalize

# define a custom kernel to compute similarity between two compounds based on their soap matrices


class SOAP_Kernel(gpytorch.kernels.kernel):
    def __init__(self) -> None:
        super().__init__()

    def forward(self, x1: np.ndarray, x2: np.ndarray, **params):
        K = REMatchKernel(
            metric="polynomial",
            degree=3,
            gamma=1,
            coef0=0,
            alpha=0.5,
            threshold=1e-6,
            normalize_kernel=True,
        )
        x1 = normalize(x1)
        x2 = normalize(x2)
        return torch.tensor(K.create(x1, x2))


class GPyTorchModel(gpytorch.models.ExactGP):
    def __init__(self, train_inputs, train_targets, likelihood):
        super().__init__(train_inputs, train_targets, likelihood)
        self.covar_module = SOAP_Kernel()
        self.mean_module = gpytorch.means.ConstantMean()

    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.multitask_multivariate_normal(mean_x, covar_x)
