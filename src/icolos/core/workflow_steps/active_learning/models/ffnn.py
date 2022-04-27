import torch
from torch import nn


class FeedForwardNet(nn.Module):
    def __init__(self) -> None:
        super(FeedForwardNet, self).__init__()
        self.net = nn.Sequential(
            nn.Linear(2048, 128),
            nn.LeakyReLU(),
            nn.Linear(128, 128),
            nn.LeakyReLU(),
            nn.Linear(128, 128),
            nn.LeakyReLU(),
            nn.Dropout(),
            nn.Linear(128, 1),
        )

    def forward(self, inputs) -> torch.Tensor:
        x = self.net(inputs)
        return x
