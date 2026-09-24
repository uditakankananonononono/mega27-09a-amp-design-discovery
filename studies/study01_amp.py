"""Study 1 - AMP-Forge: antimicrobial peptide classification benchmark + design.

Real data: reviewed UniProt KW-0929 'Antimicrobial' sequences (8-60 aa) vs
length-matched reviewed non-AMP UniProt sequences (fetched live, cached in
data/raw/). Models: PepCNN, PepGNN, and two classic baselines. Metrics:
5-fold out-of-fold AUROC/AUPRC with bootstrap 95% CIs and permutation tests.
Design: ensemble-scored simulated annealing with novelty control vs the
training corpus.
"""
from __future__ import annotations
import json
import pathlib

import numpy as np
import torch
from sklearn.metrics import roc_auc_score, average_precision_score

from pepdesign.benchmark import bootstrap_ci, permutation_test, composition_baseline, cv_evaluate
from pepdesign.encoding import combined_features
from pepdesign.models.cnn import PepCNN
from pepdesign.models.gnn import PepGNN, graph_inputs
from pepdesign.models.ensemble import train_binary, predict_scores
from pepdesign.design import anneal_design, novelty

ROOT = pathlib.Path(__file__).resolve().parent.parent
MAX_LEN = 60


def load_fasta(path):
    seqs = []
    for line in open(path):
        if line.startswith(">"):
            seqs.append("")
        elif seqs or line.strip():
            seqs[-1] += line.strip()
    return seqs


def cnn_inputs(seqs, max_len):
    return (torch.from_numpy(np.stack([combined_features(s, max_len) for s in seqs])),)


def gnn_inputs(seqs, max_len):
    return graph_inputs(seqs, max_len)


def main():
    pos = load_fasta(ROOT / "data/raw/amp_positives.fasta")
    neg = load_fasta(ROOT / "data/raw/amp_negatives.fasta")
    seqs = pos + neg
    labels = np.array([1] * len(pos) + [0] * len(neg))
    print(f"dataset: {len(pos)} AMP vs {len(neg)} non-AMP (real, UniProt reviewed)")

    results = {"dataset": {"positives": len(pos), "negatives": len(neg),
                           "source": "UniProt reviewed, keyword KW-0929, 8-60 aa"}}

    def make_tp(model_ctor, inputs_fn, epochs):
        def tp(train_seqs, train_y, test_seqs, seed):
            model = model_ctor()
            train_binary(model, inputs_fn, train_seqs, train_y, MAX_LEN,
                         epochs=epochs, seed=seed,
                         pos_weight=len(train_y) / max(2 * train_y.sum(), 1))
            return predict_scores(model, inputs_fn, test_seqs, MAX_LEN)
        return tp

    runners = {
        "PepCNN": make_tp(PepCNN, cnn_inputs, 30),
        "PepGNN": make_tp(PepGNN, gnn_inputs, 30),
        "logreg_dipeptide": lambda tr, y, te, seed: composition_baseline(tr, y, te, "logreg", seed),
        "rf_dipeptide": lambda tr, y, te, seed: composition_baseline(tr, y, te, "rf", seed),
    }
    oof = {}
    for name, tp in runners.items():
        print(f"running {name} 5-fold CV ...", flush=True)
        oof[name] = cv_evaluate(tp, seqs, labels, folds=5, seed=42)
        auc, lo, hi = bootstrap_ci(labels, oof[name], n_boot=500, seed=1)
        ap = average_precision_score(labels, oof[name])
        _, p = permutation_test(labels, oof[name], n_perm=500, seed=1)
        results[name] = {"auroc": auc, "auroc_ci95": [lo, hi], "auprc": float(ap),
                         "permutation_p": p}
        print(f"  {name}: AUROC {auc:.3f} [{lo:.3f},{hi:.3f}] AUPRC {ap:.3f} p={p:.4f}", flush=True)

    # ---- design: train final ensemble on all data, anneal with novelty control
    print("training final ensemble for design ...", flush=True)
    cnn = train_binary(PepCNN(), cnn_inputs, seqs, labels, MAX_LEN, epochs=35, seed=7,
                       pos_weight=len(labels) / (2 * labels.sum()))
    gnn = train_binary(PepGNN(), gnn_inputs, seqs, labels, MAX_LEN, epochs=35, seed=11,
                       pos_weight=len(labels) / (2 * labels.sum()))

    def ensemble_score(seq):
        s1 = predict_scores(cnn, cnn_inputs, [seq], MAX_LEN)[0]
        s2 = predict_scores(gnn, gnn_inputs, [seq], MAX_LEN)[0]
        return 0.5 * (float(s1) + float(s2))

    corpus = pos  # novelty measured against real AMPs
    designs = []
    for i in range(24):
        res = anneal_design(ensemble_score, length=25, steps=250,
                            corpus=corpus, novelty_floor=0.35, seed=100 + i)
        designs.append({"sequence": res.best_sequence, "score": res.best_score,
                        "novelty_vs_corpus": novelty(res.best_sequence, corpus)})
    designs.sort(key=lambda d: -d["score"])
    results["designs"] = designs[:10]
    for d in designs[:5]:
        print(f"  design score={d['score']:.3f} novelty={d['novelty_vs_corpus']:.2f} {d['sequence']}")

    out = ROOT / "results/study01_amp.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"wrote {out}")

    # ROC figure for the paper
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve
    plt.figure(figsize=(6, 5))
    for name in runners:
        fpr, tpr, _ = roc_curve(labels, oof[name])
        plt.plot(fpr, tpr, label=f"{name} (AUC {roc_auc_score(labels, oof[name]):.3f})")
    plt.plot([0, 1], [0, 1], "k--", lw=0.8)
    plt.xlabel("False positive rate"); plt.ylabel("True positive rate")
    plt.title("AMP-Forge: out-of-fold ROC (5-fold CV, real UniProt data)")
    plt.legend(fontsize=8); plt.tight_layout()
    plt.savefig(ROOT / "results/study01_roc.png", dpi=150)
    print("wrote results/study01_roc.png")


if __name__ == "__main__":
    main()
