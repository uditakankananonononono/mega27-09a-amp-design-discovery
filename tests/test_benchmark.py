import numpy as np
from pepdesign.benchmark import (bootstrap_ci, permutation_test,
                                 composition_baseline, cv_evaluate)


def test_bootstrap_ci_perfect_scores():
    y = [0] * 50 + [1] * 50
    s = [0.1] * 50 + [0.9] * 50
    auc, lo, hi = bootstrap_ci(y, s, n_boot=100, seed=0)
    assert auc == 1.0 and lo <= auc <= hi


def test_permutation_test_distinguishes_signal():
    rng = np.random.default_rng(0)
    y = np.array([0] * 60 + [1] * 60)
    good = y + rng.normal(0, 0.2, size=120)
    _, p_good = permutation_test(y, good, n_perm=200, seed=1)
    noise = rng.random(120)
    _, p_noise = permutation_test(y, noise, n_perm=200, seed=1)
    assert p_good < 0.05
    assert p_noise > p_good


def test_composition_baseline_shape():
    train = ["AAAA", "CCCC", "KKKK", "DDDD"]
    scores = composition_baseline(train, [0, 0, 1, 1], train, seed=0)
    assert scores.shape == (4,)
    assert ((scores >= 0) & (scores <= 1)).all()


def test_cv_evaluate_oof_alignment():
    seqs = ["AAAA"] * 10 + ["KKKK"] * 10
    labels = np.array([0] * 10 + [1] * 10)

    def tp(tr, y, te, seed):
        return np.array([0.2] * len(te))

    oof = cv_evaluate(tp, seqs, labels, folds=5, seed=0)
    assert oof.shape == (20,)
    assert np.allclose(oof, 0.2)
