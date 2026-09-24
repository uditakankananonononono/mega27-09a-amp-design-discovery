"""In-silico peptide design loop: ensemble scoring, simulated-annealing mutation,
novelty control against the training corpus, and full trajectory logging."""
from __future__ import annotations
import math
from dataclasses import dataclass, field

import numpy as np

from pepdesign.encoding import AA, clean_sequence


@dataclass
class DesignRecord:
    sequence: str
    score: float
    step: int
    accepted: bool


@dataclass
class DesignResult:
    best_sequence: str
    best_score: float
    trajectory: list = field(default_factory=list)


def kmer_set(seq: str, k: int = 3) -> set:
    s = clean_sequence(seq)
    return {s[i:i + k] for i in range(len(s) - k + 1)} if len(s) >= k else {s} if s else set()


def novelty(seq: str, corpus, k: int = 3) -> float:
    """1 - max Jaccard(k-mer set) against any corpus sequence."""
    mine = kmer_set(seq, k)
    if not mine:
        return 0.0
    best = 0.0
    for other in corpus:
        theirs = kmer_set(other, k)
        if not theirs:
            continue
        best = max(best, len(mine & theirs) / len(mine | theirs))
    return 1.0 - best


def random_peptide(length: int, rng: np.random.Generator) -> str:
    return "".join(rng.choice(list(AA), size=length))


def mutate(seq: str, rng: np.random.Generator, n_mut: int = 1) -> str:
    s = list(seq)
    for _ in range(n_mut):
        i = rng.integers(0, len(s))
        s[i] = rng.choice([a for a in AA if a != s[i]])
    return "".join(s)


def anneal_design(score_fn, length: int, steps: int = 300, seed_seq: str | None = None,
                  corpus=None, novelty_floor: float = 0.15, t0: float = 0.2,
                  cooling: float = 0.99, seed: int = 0) -> DesignResult:
    """Simulated annealing over sequence space. score_fn(seq) -> float in [0,1]."""
    rng = np.random.default_rng(seed)
    seq = seed_seq or random_peptide(length, rng)
    score = float(score_fn(seq))
    best = DesignResult(best_sequence=seq, best_score=score)
    temp = t0
    for step in range(steps):
        cand = mutate(seq, rng, n_mut=1 + int(rng.random() < 0.2))
        cscore = float(score_fn(cand))
        nov_ok = corpus is None or novelty(cand, corpus) >= novelty_floor
        accept = cscore >= score or (nov_ok and rng.random() < math.exp((cscore - score) / max(temp, 1e-6)))
        if accept and nov_ok:
            seq, score = cand, cscore
        best.trajectory.append(DesignRecord(cand, cscore, step, accept and nov_ok))
        if score > best.best_score:
            best.best_sequence, best.best_score = seq, score
        temp *= cooling
    return best
