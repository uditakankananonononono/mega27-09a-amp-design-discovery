"""Training/evaluation utilities shared by the CNN, GNN, and baselines."""
from __future__ import annotations
import numpy as np
import torch
import torch.nn as nn


def train_binary(model: nn.Module, make_inputs, seqs, labels, max_len: int,
                 epochs: int = 40, lr: float = 1e-3, batch_size: int = 64,
                 pos_weight: float = 1.0, seed: int = 0, verbose: bool = False):
    torch.manual_seed(seed)
    x_all = make_inputs(seqs, max_len)
    y = torch.tensor(np.asarray(labels, dtype=np.float32))
    if isinstance(x_all, tuple):
        tensors = x_all
        n = y.shape[0]
        idx = torch.randperm(n)
        opt = torch.optim.Adam(model.parameters(), lr=lr)
        lossf = nn.BCEWithLogitsLoss(pos_weight=torch.tensor(pos_weight))
        model.train()
        for ep in range(epochs):
            for s in range(0, n, batch_size):
                b = idx[s:s + batch_size]
                logits = model(*(t[b] for t in tensors))
                loss = lossf(logits, y[b])
                opt.zero_grad(); loss.backward(); opt.step()
            if verbose and ep % 10 == 0:
                print(f"  epoch {ep}: loss={loss.item():.4f}")
    else:
        raise TypeError("make_inputs must return a tuple of tensors")
    model.eval()
    return model


@torch.no_grad()
def predict_scores(model: nn.Module, make_inputs, seqs, max_len: int,
                   batch_size: int = 256) -> np.ndarray:
    x_all = make_inputs(seqs, max_len)
    out = []
    for s in range(0, len(seqs), batch_size):
        logits = model(*(t[s:s + batch_size] for t in x_all))
        out.append(torch.sigmoid(logits).numpy())
    return np.concatenate(out)
