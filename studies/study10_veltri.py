"""Study 10 - head-to-head on the AMP Scanner Vr.2 benchmark (Veltri et al. 2018).

Data: original-dataset from https://github.com/dan-veltri/amp-scanner-v2
(AMP.tr/eval/te + DECOY.tr/eval/te, the exact splits distributed with the
published model). Same splits, our models trained from scratch: PepCNN,
PepGNN (shared-adj), logreg/RF dipeptide, CNN+RF ensemble. Metrics on the
published test split (712 AMP / 712 DECOY): accuracy, sensitivity,
specificity, MCC, AUROC - compared against the published AMP Scanner Vr.2
numbers on this benchmark (see paper, amp scanner vr2 news page: 10-fold CV
~90.4% acc / 0.81 MCC / 0.962 auROC; holdout test ~0.91 acc).
"""
from __future__ import annotations
import json, pathlib

import numpy as np
from sklearn.metrics import (roc_auc_score, average_precision_score,
                             accuracy_score, matthews_corrcoef)

from pepdesign.benchmark import bootstrap_ci
from pepdesign.encoding import combined_features
from pepdesign.models.cnn import PepCNN
from pepdesign.models.gnn import PepGNN
from pepdesign.models.ensemble import (train_binary, predict_scores,
                                       train_binary_shared_adj,
                                       predict_scores_shared_adj)
from pepdesign.benchmark import composition_baseline

ROOT = pathlib.Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw/veltri"
MAX_LEN = 150
CNN_EPOCHS = 12


def load_fasta(path):
    seqs, cur = [], None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if cur: seqs.append(cur)
                cur = ""
            elif line:
                cur = (cur or "") + line
    if cur: seqs.append(cur)
    return seqs


def main():
    tr_s = load_fasta(RAW / "AMP.tr.fa") + load_fasta(RAW / "DECOY.tr.fa")
    tr_y = np.array([1] * 712 + [0] * 712)
    va_s = load_fasta(RAW / "AMP.eval.fa") + load_fasta(RAW / "DECOY.eval.fa")
    va_y = np.array([1] * 354 + [0] * 354)
    te_s = load_fasta(RAW / "AMP.te.fa") + load_fasta(RAW / "DECOY.te.fa")
    te_y = np.array([1] * 712 + [0] * 712)
    print(f"splits: train {len(tr_s)} eval {len(va_s)} test {len(te_s)}", flush=True)
    print(f"len range: train {min(map(len,tr_s))}-{max(map(len,tr_s))} test {min(map(len,te_s))}-{max(map(len,te_s))}", flush=True)

    results = {"benchmark": "AMP Scanner Vr.2 original-dataset (Veltri et al. 2018)",
               "splits": {"train": len(tr_s), "eval": len(va_s), "test": len(te_s)},
               "max_len": MAX_LEN,
               "published_amp_scanner_vr2": {
                   "note": "10-fold CV on this dataset per dveltri.com/ascan/v2 news",
                   "accuracy": 0.904, "mcc": 0.812, "auroc": 0.962}}

    import torch
    def cnn_inputs(seqs, max_len):
        return (torch.stack([torch.from_numpy(combined_features(x, max_len)) for x in seqs]),)
    def report(name, scores):
        auc = roc_auc_score(te_y, scores)
        _, lo, hi = bootstrap_ci(te_y, scores, n_boot=500, seed=1)
        ap = average_precision_score(te_y, scores)
        pred = (scores >= 0.5).astype(int)
        acc = accuracy_score(te_y, pred)
        mcc = matthews_corrcoef(te_y, pred)
        tp = int(((pred == 1) & (te_y == 1)).sum()); fn = int(((pred == 0) & (te_y == 1)).sum())
        tn = int(((pred == 0) & (te_y == 0)).sum()); fp = int(((pred == 1) & (te_y == 0)).sum())
        sens = tp / (tp + fn); spec = tn / (tn + fp)
        results[name] = {"accuracy": acc, "sensitivity": sens, "specificity": spec,
                         "mcc": float(mcc), "auroc": auc, "auroc_ci95": [lo, hi],
                         "auprc": float(ap)}
        print(f"  {name}: acc {acc:.3f} sens {sens:.3f} spec {spec:.3f} MCC {mcc:.3f} "
              f"AUROC {auc:.3f} [{lo:.3f},{hi:.3f}] AUPRC {ap:.3f}", flush=True)
        return scores

    # Final protocol: train on tr, early validation signal from eval, test on te.
    pw = len(tr_y) / max(2 * tr_y.sum(), 1)
    print("training PepCNN ...", flush=True)
    cnn = train_binary(PepCNN(), cnn_inputs, tr_s, tr_y, MAX_LEN,
                       epochs=CNN_EPOCHS, batch_size=64, pos_weight=pw, seed=7, verbose=True)
    s_cnn = report("PepCNN", predict_scores(cnn, cnn_inputs, te_s, MAX_LEN))

    print("training PepGNN (shared adj) ...", flush=True)
    gnn = train_binary_shared_adj(PepGNN(), tr_s, tr_y, MAX_LEN, epochs=10,
                                  batch_size=64, pos_weight=pw, seed=11, verbose=True)
    s_gnn = report("PepGNN", predict_scores_shared_adj(gnn, te_s, MAX_LEN))

    print("training rf_dipeptide ...", flush=True)
    s_rf = report("rf_dipeptide", composition_baseline(tr_s, tr_y, te_s, "rf", 42))
    print("training logreg_dipeptide ...", flush=True)
    report("logreg_dipeptide", composition_baseline(tr_s, tr_y, te_s, "logreg", 42))

    report("ensemble_cnn_rf", 0.5 * (np.asarray(s_cnn) + np.asarray(s_rf)))

    out = ROOT / "results/study10_veltri.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
