import gpytorch
from dscribe.kernels import REMatchKernel
import torch
import numpy as np
from sklearn.preprocessing import normalize

# define a custom kernel to compute similarity between two compounds based on their soap matrices


class SOAP_Kernel(gpytorch.kernels.Kernel):
    def __init__(
        self,
    ) -> None:
        super().__init__(active_dims=[0])

    def forward(self, x1: torch.Tensor, x2: torch.Tensor = None):

        soap_k = REMatchKernel(
            metric="polynomial",
            degree=3,
            gamma=1,
            coef0=0,
            alpha=0.5,
            threshold=1e-6,
            normalize_kernel=True,
        )
        return soap_k.create(x1, x2)


class SOAP_GP(gpytorch.models.ExactGP):
    def __init__(self, likelihood, noise_init=None):
        # detail: We don't set train_inputs and train_targets here skorch because
        # will take care of that.
        super().__init__(train_inputs=None, train_targets=None, likelihood=likelihood)
        self.mean_module = gpytorch.means.ConstantMean()
        self.covar_module = SOAP_Kernel()

    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)
