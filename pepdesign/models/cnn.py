"""PepCNN: 1D convolutional peptide classifier over one-hot+physicochemical channels."""
from __future__ import annotations
import torch
import torch.nn as nn


class PepCNN(nn.Module):
    def __init__(self, in_channels: int = 24, num_classes: int = 1,
                 hidden: int = 32, kernel_sizes=(3, 5, 7), dropout: float = 0.3):
        super().__init__()
        self.convs = nn.ModuleList([
            nn.Sequential(
                nn.Conv1d(in_channels, hidden, k, padding=k // 2),
                nn.ReLU(),
                nn.Conv1d(hidden, hidden, k, padding=k // 2),
                nn.ReLU(),
            ) for k in kernel_sizes])
        self.head = nn.Sequential(
            nn.Linear(hidden * len(kernel_sizes), 64),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(64, num_classes),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        # x: (batch, max_len, in_channels) -> conv over length axis
        x = x.transpose(1, 2)
        pooled = [torch.amax(c(x), dim=2) for c in self.convs]  # global max pool per branch
        return self.head(torch.cat(pooled, dim=1)).squeeze(-1)


def count_parameters(model: nn.Module) -> int:
    return sum(p.numel() for p in model.parameters() if p.requires_grad)
