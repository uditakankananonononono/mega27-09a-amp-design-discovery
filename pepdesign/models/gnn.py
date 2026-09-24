"""PepGNN: pure-PyTorch graph convolutional peptide classifier over the residue
chain graph (sequential window adjacency + self loops)."""
from __future__ import annotations
import numpy as np
import torch
import torch.nn as nn


class GraphConv(nn.Module):
    """GCN layer: H' = activation(A_norm H W + b). A_norm precomputed per sample."""

    def __init__(self, in_dim: int, out_dim: int):
        super().__init__()
        self.weight = nn.Parameter(torch.empty(in_dim, out_dim))
        nn.init.xavier_uniform_(self.weight)
        self.bias = nn.Parameter(torch.zeros(out_dim))

    def forward(self, h: torch.Tensor, adj: torch.Tensor) -> torch.Tensor:
        # h: (batch, n, in_dim); adj: (batch, n, n)
        return torch.relu(torch.matmul(adj, h) @ self.weight + self.bias)


class PepGNN(nn.Module):
    def __init__(self, node_dim: int = 4, hidden: int = 32, layers: int = 3,
                 num_classes: int = 1, dropout: float = 0.3):
        super().__init__()
        dims = [node_dim] + [hidden] * layers
        self.convs = nn.ModuleList([GraphConv(dims[i], dims[i + 1]) for i in range(layers)])
        self.head = nn.Sequential(
            nn.Linear(hidden, 64), nn.ReLU(), nn.Dropout(dropout), nn.Linear(64, num_classes))

    def forward(self, node_feats: torch.Tensor, adj: torch.Tensor) -> torch.Tensor:
        h = node_feats
        for conv in self.convs:
            h = conv(h, adj)
        pooled = h.mean(dim=1)  # mean pool over residues
        return self.head(pooled).squeeze(-1)


def graph_inputs(seqs, max_len: int) -> tuple[torch.Tensor, torch.Tensor]:
    from pepdesign.encoding import physchem, chain_adjacency
    feats = np.stack([physchem(s, max_len) for s in seqs])
    adjs = np.stack([chain_adjacency(max_len) for _ in seqs])
    return torch.from_numpy(feats), torch.from_numpy(adjs)


def graph_inputs_shared(seqs, max_len: int) -> tuple[torch.Tensor, torch.Tensor]:
    """Memory-efficient variant: node feats stacked per sample, ONE shared
    chain adjacency (max_len, max_len) broadcast by matmul in GraphConv.
    Per-sample adjacency copies cost O(n * max_len^2); at 150 aa and 50k
    sequences that is ~4.5 GB, which does not fit this sandbox."""
    from pepdesign.encoding import physchem, chain_adjacency
    feats = np.stack([physchem(s, max_len) for s in seqs])
    return torch.from_numpy(feats), torch.from_numpy(chain_adjacency(max_len))
