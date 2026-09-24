"""Sequence encodings for peptides: one-hot, physicochemical, k-mer composition,
and residue-graph construction for the GNN."""
from __future__ import annotations
import numpy as np

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {a: i for i, a in enumerate(AA)}

# Kyte-Doolittle hydrophobicity
HYDROPHOBICITY = dict(zip(AA, [1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8,
                               1.9, -3.5, -1.6, -3.5, -4.5, -0.8, -0.7, 4.2, -0.9, -1.3]))
# Net charge at physiological pH
CHARGE = {a: 0.0 for a in AA}
CHARGE.update({"K": 1.0, "R": 1.0, "H": 0.1, "D": -1.0, "E": -1.0})
# Approximate residue volume (A^3)
VOLUME = dict(zip(AA, [88.6, 108.5, 111.1, 138.4, 189.9, 60.1, 153.2, 166.7, 168.6, 166.7,
                       166.7, 124.9, 135.8, 175.4, 112.7, 89.0, 116.1, 227.8, 193.6, 140.0]))

PHYSCHEM_DIM = 4  # hydrophobicity, charge, volume, is_aromatic
AROMATIC = {"F", "W", "Y"}


def clean_sequence(seq: str) -> str:
    return "".join(a for a in seq.upper() if a in AA_INDEX)


def one_hot(seq: str, max_len: int) -> np.ndarray:
    """(max_len, 20) one-hot matrix, zero-padded."""
    out = np.zeros((max_len, len(AA)), dtype=np.float32)
    for i, a in enumerate(clean_sequence(seq)[:max_len]):
        out[i, AA_INDEX[a]] = 1.0
    return out


def physchem(seq: str, max_len: int) -> np.ndarray:
    """(max_len, PHYSCHEM_DIM) per-residue physicochemical matrix."""
    out = np.zeros((max_len, PHYSCHEM_DIM), dtype=np.float32)
    for i, a in enumerate(clean_sequence(seq)[:max_len]):
        out[i] = [HYDROPHOBICITY[a] / 4.5, CHARGE[a], VOLUME[a] / 230.0, 1.0 if a in AROMATIC else 0.0]
    return out


def combined_features(seq: str, max_len: int) -> np.ndarray:
    """(max_len, 24): one-hot concatenated with physicochemical features."""
    return np.concatenate([one_hot(seq, max_len), physchem(seq, max_len)], axis=1)


def composition_vector(seq: str, k: int = 2) -> np.ndarray:
    """Normalized k-mer composition vector over the canonical alphabet."""
    s = clean_sequence(seq)
    dim = len(AA) ** k
    vec = np.zeros(dim, dtype=np.float64)
    if len(s) < k:
        return vec
    for i in range(len(s) - k + 1):
        idx = 0
        for a in s[i:i + k]:
            idx = idx * len(AA) + AA_INDEX[a]
        vec[idx] += 1.0
    total = vec.sum()
    return vec / total if total > 0 else vec


def chain_adjacency(length: int, window: int = 3) -> np.ndarray:
    """Symmetric normalized adjacency for a residue chain graph with a local
    sequence window (|i-j| <= window), self-loops included, D^-1/2 (A+I) D^-1/2."""
    n = max(length, 1)
    a = np.zeros((n, n), dtype=np.float32)
    for i in range(n):
        for j in range(max(0, i - window), min(n, i + window + 1)):
            a[i, j] = 1.0
    a += np.eye(n, dtype=np.float32)
    d = a.sum(axis=1)
    d[d == 0] = 1.0
    d_inv_sqrt = 1.0 / np.sqrt(d)
    return (a * d_inv_sqrt[:, None]) * d_inv_sqrt[None, :]
