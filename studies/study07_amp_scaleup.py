"""Study 7 - scale-up AMP benchmark on the homology-aware split.

Data: data/processed_large/{train,val,test}.csv from build_large_dataset.py
(27k reviewed UniProt KW-0929 AMPs vs length-matched reviewed non-AMP,
cluster-aware 80/10/10 split over shared-10-mer union-find clusters, so no
near-duplicate crosses splits). Models: PepCNN, PepGNN (shared-adjacency
path, memory-bounded), logreg/RF dipeptide baselines. Metrics on the held-out
test split: AUROC with bootstrap 95% CI, AUPRC, permutation p-value.
Design: ensemble-trained annealing with novelty control at MAX_LEN=150.
"""
from __future__ import annotations
import csv, json, pathlib

import numpy as np
import torch
from sklearn.metrics import roc_auc_score, average_precision_score

from pepdesign.benchmark import bootstrap_ci, permutation_test, composition_baseline
from pepdesign.encoding import combined_features
from pepdesign.models.cnn import PepCNN
from pepdesign.models.gnn import PepGNN
from pepdesign.models.ensemble import (train_binary, predict_scores,
                                       train_binary_shared_adj,
                                       predict_scores_shared_adj)
from pepdesign.design import anneal_design, novelty

ROOT = pathlib.Path(__file__).resolve().parent.parent
MAX_LEN = 150
TRAIN_CAP = 30000          # balanced cap; fits RAM (30k x 150 x 24 floats ~ 430 MB)
CNN_EPOCHS = 8
GNN_EPOCHS = 8


def load_split(name):
    seqs, labels = [], []
    with open(ROOT / "data/processed_large" / f"{name}.csv") as f:
        for row in csv.DictReader(f):
            seqs.append(row["sequence"]); labels.append(int(row["label"]))
    return seqs, np.array(labels)


def cnn_inputs(seqs, max_len):
    return (torch.from_numpy(np.stack([combined_features(s, max_len) for s in seqs])),)


def main():
    tr_s, tr_y = load_split("train")
    va_s, va_y = load_split("val")
    te_s, te_y = load_split("test")
    print(f"split sizes: train {len(tr_s)} val {len(va_s)} test {len(te_s)}", flush=True)

    # Balanced subsample of train for the CNN/GNN (baselines use vectors, keep full train
    # only if cheap - dipeptide logreg/RF on 40k is fine too, so subsample all equally).
    rng = np.random.default_rng(27)
    pos_idx = np.where(tr_y == 1)[0]; neg_idx = np.where(tr_y == 0)[0]
    half = min(TRAIN_CAP // 2, len(pos_idx), len(neg_idx))
    keep = np.concatenate([rng.choice(pos_idx, half, replace=False),
                           rng.choice(neg_idx, half, replace=False)])
    rng.shuffle(keep)
    sub_s = [tr_s[i] for i in keep]; sub_y = tr_y[keep]
    print(f"training subsample: {len(sub_s)} ({int(sub_y.sum())} pos)", flush=True)

    pw = len(sub_y) / max(2 * sub_y.sum(), 1)
    results = {"split": {"train": len(tr_s), "val": len(va_s), "test": len(te_s),
                         "train_subsample": len(sub_s)},
               "max_len": MAX_LEN,
               "source": "UniProt reviewed KW-0929 (27k AMP) vs reviewed non-AMP, "
                         "cluster-aware 10-mer split"}

    def report(name, scores):
        auc, lo, hi = bootstrap_ci(te_y, scores, n_boot=500, seed=1)
        ap = average_precision_score(te_y, scores)
        _, p = permutation_test(te_y, scores, n_perm=500, seed=1)
        results[name] = {"auroc": auc, "auroc_ci95": [lo, hi], "auprc": float(ap),
                         "permutation_p": p}
        print(f"  {name}: AUROC {auc:.3f} [{lo:.3f},{hi:.3f}] AUPRC {ap:.3f} p={p:.4f}",
              flush=True)

    print("training PepCNN ...", flush=True)
    cnn = train_binary(PepCNN(), cnn_inputs, sub_s, sub_y, MAX_LEN,
                       epochs=CNN_EPOCHS, batch_size=128, pos_weight=pw, seed=7, verbose=True)
    report("PepCNN", predict_scores(cnn, cnn_inputs, te_s, MAX_LEN))

    print("training PepGNN (shared adj) ...", flush=True)
    gnn = train_binary_shared_adj(PepGNN(), sub_s, sub_y, MAX_LEN,
                                  epochs=GNN_EPOCHS, batch_size=128, pos_weight=pw,
                                  seed=11, verbose=True)
    report("PepGNN", predict_scores_shared_adj(gnn, te_s, MAX_LEN))

    for name, kind in [("logreg_dipeptide", "logreg"), ("rf_dipeptide", "rf")]:
        print(f"training {name} ...", flush=True)
        report(name, composition_baseline(sub_s, sub_y, te_s, kind, 42))

    # ---- design at scale: anneal with the trained ensemble, novelty vs train positives
    print("design screen (ensemble annealing) ...", flush=True)
    def ensemble_score(seq):
        s1 = predict_scores(cnn, cnn_inputs, [seq], MAX_LEN)[0]
        s2 = predict_scores_shared_adj(gnn, [seq], MAX_LEN)[0]
        return 0.5 * (float(s1) + float(s2))

    corpus = [tr_s[i] for i in pos_idx[:3000]]  # novelty reference subsample
    designs = []
    for i in range(24):
        res = anneal_design(ensemble_score, length=25, steps=250,
                            corpus=corpus, novelty_floor=0.35, seed=500 + i)
        designs.append({"sequence": res.best_sequence, "score": res.best_score,
                        "novelty_vs_corpus": novelty(res.best_sequence, corpus)})
        print(f"  design {i}: score={res.best_score:.3f}", flush=True)
    designs.sort(key=lambda d: -d["score"])
    results["designs"] = designs[:10]

    torch.save({"state_dict": cnn.state_dict(), "max_len": MAX_LEN},
               ROOT / "results/study07_cnn.pt")
    torch.save({"state_dict": gnn.state_dict(), "max_len": MAX_LEN},
               ROOT / "results/study07_gnn.pt")
    out = ROOT / "results/study07_scaleup.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
