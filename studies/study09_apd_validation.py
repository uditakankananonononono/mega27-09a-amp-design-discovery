"""Study 9 - external validation on APD6 (2024 natural AMPs).

Scores the APD6 2024 natural AMP list (3,306 entries, aps.unmc.edu
downloads, APD3 successor) with the study07 PepCNN + PepGNN checkpoints.
APD entries sharing an exact 10-mer with any large-split training
sequence are excluded (decontamination). Negatives are the held-out
test-split non-AMPs. Reports CNN / GNN / ensemble AUROC + AUPRC with
bootstrap 95% CIs. APD6 is external to the UniProt KW-0929 training
source, so this measures cross-database generalization.
"""
from __future__ import annotations
import csv, json, pathlib

import numpy as np
import torch
from sklearn.metrics import roc_auc_score, average_precision_score

from pepdesign.encoding import combined_features
from pepdesign.models.cnn import PepCNN
from pepdesign.models.gnn import PepGNN
from pepdesign.models.ensemble import predict_scores, predict_scores_shared_adj

ROOT = pathlib.Path(__file__).resolve().parent.parent
KMER = 10
MAX_LEN_CAP = 150
BOOT = 2000
SEED = 27
STANDARD = set("ACDEFGHIKLMNPQRSTVWY")


def parse_fasta(path):
    recs, header, seq = [], None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if header is not None:
                recs.append((header, "".join(seq)))
            header, seq = line[1:], []
        elif line:
            seq.append(line)
    if header is not None:
        recs.append((header, "".join(seq)))
    return recs


def cnn_inputs(seqs, max_len):
    return (torch.from_numpy(np.stack([combined_features(s, max_len) for s in seqs])),)


def kmers(seq, k):
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def bootstrap_ci(y, s, fn, rng, n=BOOT):
    vals = []
    y = np.asarray(y); s = np.asarray(s)
    for _ in range(n):
        idx = rng.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) < 2:
            continue
        vals.append(fn(y[idx], s[idx]))
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def main():
    rng = np.random.default_rng(SEED)
    cnn_ckpt = torch.load(ROOT / "results/study07_cnn.pt", weights_only=False)
    gnn_ckpt = torch.load(ROOT / "results/study07_gnn.pt", weights_only=False)
    max_len = cnn_ckpt["max_len"]
    cnn = PepCNN(); cnn.load_state_dict(cnn_ckpt["state_dict"]); cnn.eval()
    gnn = PepGNN(); gnn.load_state_dict(gnn_ckpt["state_dict"]); gnn.eval()

    # training 10-mer index for decontamination
    train_kmers = set()
    with open(ROOT / "data/processed_large/train.csv") as f:
        for row in csv.DictReader(f):
            train_kmers |= kmers(row["sequence"], KMER)
    print(f"train {KMER}-mers indexed: {len(train_kmers)}", flush=True)

    # APD6 positives: standard residues, 10-150 aa, dedup, drop the note header
    recs = parse_fasta(ROOT / "data/raw/apd6_natural_2024.fasta")
    seen, pos = set(), []
    for h, s in recs:
        s = s.upper()
        if not (10 <= len(s) <= MAX_LEN_CAP) or not set(s) <= STANDARD or s in seen:
            continue
        seen.add(s)
        pos.append((h, s))
    n_raw = len(pos)
    pos = [(h, s) for h, s in pos if not (kmers(s, KMER) & train_kmers)]
    print(f"APD6 usable: {n_raw}, after decontamination: {len(pos)}", flush=True)

    # negatives: held-out test-split non-AMPs
    neg = []
    with open(ROOT / "data/processed_large/test.csv") as f:
        for row in csv.DictReader(f):
            if row["label"] == "0":
                neg.append(row["sequence"])
    print(f"test-split non-AMP negatives: {len(neg)}", flush=True)

    seqs = [s for _, s in pos] + neg
    y = np.array([1] * len(pos) + [0] * len(neg))
    s_cnn = predict_scores(cnn, cnn_inputs, seqs, max_len)
    s_gnn = predict_scores_shared_adj(gnn, seqs, max_len)
    s_ens = 0.5 * (s_cnn + s_gnn)

    out = {"apd6_positives_raw": n_raw, "apd6_positives_decontaminated": len(pos),
           "negatives_test_split": len(neg), "models": {}}
    for name, s in [("PepCNN", s_cnn), ("PepGNN", s_gnn), ("ensemble", s_ens)]:
        auroc = float(roc_auc_score(y, s))
        auprc = float(average_precision_score(y, s))
        lo, hi = bootstrap_ci(y, s, roc_auc_score, rng)
        out["models"][name] = {"auroc": auroc, "auroc_ci95": [lo, hi], "auprc": auprc}
        print(f"{name}: AUROC {auroc:.4f} [{lo:.4f}, {hi:.4f}] AUPRC {auprc:.4f}", flush=True)

    # sensitivity of APD6 hits: fraction scoring above the test-split 95th pct
    test_scores = []
    res_path = ROOT / "results/study07_amp_scaleup.json"
    if res_path.exists():
        out["note"] = "sensitivity uses ensemble score distribution on APD6"
    thr = float(np.quantile(s_ens[len(pos):], 0.95))
    out["apd6_sensitivity_at_test_fpr5"] = float(np.mean(s_ens[:len(pos)] > thr))
    print(f"APD6 sensitivity @ 5% FPR: {out['apd6_sensitivity_at_test_fpr5']:.4f}", flush=True)

    path = ROOT / "results/study09_apd_validation.json"
    path.write_text(json.dumps(out, indent=2))
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
