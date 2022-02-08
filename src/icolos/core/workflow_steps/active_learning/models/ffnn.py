try:
    import torch
    from torch import nn
except ImportError:
    print(
        "Warning: PyTorch imports failed - ensure you have pytorch + pytorch lightning installed in your environment if you are expecting to use them!"
    )

"""
Implement a FFNN in PyTorch Lightning
"""


class FeedForwardNet(nn.Module):
    def __init__(self) -> None:
        super(FeedForwardNet, self).__init__()
        self.net = nn.Sequential(
            nn.Linear(2048, 256),
            nn.LeakyReLU(),
            nn.Linear(256, 1),
        )

    def forward(self, inputs) -> torch.Tensor:
        x = self.net(inputs)
        return x
