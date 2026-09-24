import numpy as np
from pepdesign.design import (anneal_design, novelty, random_peptide, mutate,
                              kmer_set)
from pepdesign.data.fixtures import AMP_FIXTURE


def score_fn(seq):
    # synthetic: reward lysine content
    return seq.count("K") / len(seq)


def test_anneal_improves_score():
    res = anneal_design(score_fn, length=12, steps=300, seed=3)
    # random expectation is ~0.05 (1/20 residues are K); annealing must beat it clearly
    assert res.best_score > 0.2
    assert res.best_score >= score_fn(res.trajectory[0].sequence)
    assert len(res.trajectory) == 300


def test_novelty_self_vs_distinct():
    assert novelty(AMP_FIXTURE[0], [AMP_FIXTURE[0]]) == 0.0
    assert novelty("KKKKKKKK", AMP_FIXTURE) > 0.5


def test_kmer_set_and_mutate():
    ks = kmer_set("ACDEFG", 3)
    assert "ACD" in ks and "EFG" in ks
    rng = np.random.default_rng(0)
    m = mutate("AAAAAAAAAA", rng, n_mut=2)
    assert len(m) == 10 and m != "AAAAAAAAAA"


def test_novelty_floor_blocks_copies():
    res = anneal_design(score_fn, length=10, steps=50,
                        corpus=AMP_FIXTURE, novelty_floor=0.95, seed=5)
    assert novelty(res.best_sequence, AMP_FIXTURE) >= 0.95 or res.best_score >= score_fn(res.trajectory[0].sequence)
