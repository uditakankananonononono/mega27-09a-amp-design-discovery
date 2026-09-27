"""Study 22 - Q1 reviewed-only clean rebuild (addendum-4, locked 2026-09-27;
amendment A4.1 disclosed before use). Same study19 ensemble (CNN7 .20/CNN27 .20/
GNN .15/RF .20/HGB .25, threshold 0.5) and rng(11) random 10-fold protocol,
retrained on the CLEAN dataset: 2,973 unique reviewed (Swiss-Prot) KW-0929
positives + 2,973 matched reviewed negatives (species/genus/global-length
fallback per A4.1). Locked reporting: metrics vs the confounded scale-up
(AUROC 0.921 / APD6 0.926); the delta IS the headline; no success threshold;
any direction reported as-is.

  PYTHONPATH=. python3 studies/study22_q1_reviewed_rebuild.py --folds 1-10
  PYTHONPATH=. python3 studies/study22_q1_reviewed_rebuild.py --finalize
"""
from __future__ import annotations
import argparse, json, math, pathlib, statistics as st
import numpy as np
from studies.study14_feb2020_full import load_fasta_seqs
from studies.study19_feb2020_v2 import run_fold

ROOT = pathlib.Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw"
PART = ROOT / "results/study22_q1_partial"
FINAL = ROOT / "results/study22_q1_reviewed_rebuild.json"
FOLDS = 10

def load_data():
    pos = load_fasta_seqs(RAW / "q1_reviewed_positives.fasta")
    neg = load_fasta_seqs(RAW / "q1_matched_negatives.fasta")
    assert len(pos) == len(neg), (len(pos), len(neg))
    seqs = pos + neg
    y = np.array([1]*len(pos) + [0]*len(neg))
    return seqs, y

def make_splits(n):
    rng = np.random.RandomState(11)
    perm = rng.permutation(n)
    return [(np.concatenate([perm[:i*n//FOLDS], perm[(i+1)*n//FOLDS:]]),
             perm[i*n//FOLDS:(i+1)*n//FOLDS]) for i in range(FOLDS)]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", type=str, default="1-10")
    ap.add_argument("--finalize", action="store_true")
    args = ap.parse_args()
    seqs, y = load_data()
    if args.finalize:
        rows = [json.load(open(p)) for p in sorted(PART.glob("fold_*.json"),
                key=lambda p: int(p.name.split('_')[1].split('.')[0]))]
        assert len(rows) == FOLDS, len(rows)
        out = {"study": "study22_q1_reviewed_rebuild", "n_pos": int(y.sum()),
               "n_neg": int((1-y).sum()), "protocol": "study19 ensemble, rng(11) random 10-fold",
               "reference_confounded": {"scaleup_auroc": 0.921, "apd6_auroc": 0.926},
               "folds": rows}
        for k in ["acc", "sens", "spec", "mcc", "auc"]:
            v = [r[k] for r in rows]; m, s = st.mean(v), st.stdev(v)
            out[k] = {"folds": [round(x,4) for x in v], "mean": round(m,4), "sd": round(s,4),
                      "lower95_one_sided": round(m - 1.833*s/math.sqrt(FOLDS), 4)}
        out["headline_delta_vs_confounded_scaleup_auroc"] = round(out["auc"]["mean"] - 0.921, 4)
        json.dump(out, open(FINAL, "w"), indent=1)
        print("FINAL:", json.dumps({k: out[k]["mean"] for k in ["acc","mcc","auc"]}),
              "delta_vs_0.921:", out["headline_delta_vs_confounded_scaleup_auroc"])
        return
    lo, hi = (int(x) for x in args.folds.split("-"))
    PART.mkdir(exist_ok=True)
    splits = make_splits(len(seqs))
    for f in range(lo, hi+1):
        pf = PART / f"fold_{f}.json"
        if pf.exists():
            print(f"fold {f}: exists, skip", flush=True); continue
        row = run_fold(f, seqs, y, splits)
        json.dump(row, open(pf, "w"))
        print(f"fold {f}: acc {row['acc']:.4f} mcc {row['mcc']:.4f} auc {row['auc']:.4f} ({row['seconds']:.1f}s)", flush=True)

if __name__ == "__main__":
    main()
