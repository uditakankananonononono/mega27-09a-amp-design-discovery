"""Benchmark harness: cross-validation, AUROC/AUPRC with bootstrap CIs,
permutation tests, and classic baselines (composition logistic regression,
random forest)."""
from __future__ import annotations
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import StratifiedKFold

from pepdesign.encoding import composition_vector


def bootstrap_ci(y_true, scores, metric=roc_auc_score, n_boot: int = 500,
                 alpha: float = 0.05, seed: int = 0):
    rng = np.random.default_rng(seed)
    y_true = np.asarray(y_true); scores = np.asarray(scores)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(y_true), len(y_true))
        if len(np.unique(y_true[idx])) < 2:
            continue
        vals.append(metric(y_true[idx], scores[idx]))
    lo, hi = np.percentile(vals, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return float(metric(y_true, scores)), float(lo), float(hi)


def permutation_test(y_true, scores, metric=roc_auc_score, n_perm: int = 500, seed: int = 0):
    """One-sided permutation p-value for the observed metric under label shuffle."""
    rng = np.random.default_rng(seed)
    y_true = np.asarray(y_true); scores = np.asarray(scores)
    obs = metric(y_true, scores)
    count = 1  # continuity correction
    for _ in range(n_perm):
        perm = rng.permutation(y_true)
        if metric(perm, scores) >= obs:
            count += 1
    return obs, count / (n_perm + 1)


def composition_baseline(train_seqs, train_y, test_seqs, method: str = "logreg", seed: int = 0):
    xtr = np.stack([composition_vector(s) for s in train_seqs])
    xte = np.stack([composition_vector(s) for s in test_seqs])
    if method == "logreg":
        clf = LogisticRegression(max_iter=500, random_state=seed)
    elif method == "rf":
        clf = RandomForestClassifier(n_estimators=100, random_state=seed, n_jobs=1)
    else:
        raise ValueError(f"unknown baseline method: {method}")
    clf.fit(xtr, train_y)
    return clf.predict_proba(xte)[:, 1]


def cv_evaluate(train_predict, seqs, labels, folds: int = 5, seed: int = 0):
    """train_predict(train_seqs, train_y, test_seqs, fold_seed) -> test scores.
    Returns out-of-fold scores aligned with `seqs`."""
    skf = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
    labels = np.asarray(labels)
    oof = np.zeros(len(seqs), dtype=np.float64)
    for k, (tr, te) in enumerate(skf.split(seqs, labels)):
        train_seqs = [seqs[i] for i in tr]
        test_seqs = [seqs[i] for i in te]
        oof[te] = train_predict(train_seqs, labels[tr], test_seqs, seed + k)
    return oof
