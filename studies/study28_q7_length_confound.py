"""Study 28 - Q7 length-confound test (addendum 7, locked 2026-09-27 BEFORE
scoring). Exact-length 1:1 matched reviewed-only subset of the study22 locked
Q1 dataset, retrained with the identical study19 ensemble and rng(11) random
10-fold protocol. The delta vs study22 (AUROC 0.9633 / MCC 0.803) IS the
length-confound estimate. No success threshold; a drop is the expected,
citable outcome; any direction reported as-is. Pre-locked fallback: if exact
matching yields < 1,000 pairs, widen to length +/- 1 aa (disclosed in output).
Falsifier: label-shuffle arm (rng(11)) must collapse to chance (|MCC| < 0.1).

  PYTHONPATH=. python3 studies/study28_q7_length_confound.py --folds 1-10
  PYTHONPATH=. python3 studies/study28_q7_length_confound.py --falsifier --folds 1-10
  PYTHONPATH=. python3 studies/study28_q7_length_confound.py --finalize
"""
from __future__ import annotations
import argparse, json, math, pathlib, statistics as st
from collections import defaultdict
import numpy as np
from studies.study14_feb2020_full import load_fasta_seqs
from studies.study19_feb2020_v2 import run_fold

ROOT = pathlib.Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw"
PART = ROOT / "results/study28_q7_partial"
FPART = ROOT / "results/study28_q7_falsifier_partial"
FINAL = ROOT / "results/study28_q7_length_confound.json"
FOLDS = 10
MIN_PAIRS = 1000

def match_exact(pos, neg, widen):
    rng = np.random.RandomState(11)
    porder = rng.permutation(len(pos))
    by_len = defaultdict(list)
    for j, s in enumerate(neg):
        by_len[len(s)].append(j)
    for L in by_len:
        by_len[L] = list(rng.permutation(by_len[L]))
    pairs, dropped = [], 0
    for i in porder:
        L = len(pos[i])
        cand = by_len.get(L, [])
        if not cand and widen:
            cand = by_len.get(L-1, []) or by_len.get(L+1, [])
        if cand:
            pairs.append((i, cand.pop()))
        else:
            dropped += 1
    return pairs, dropped

def load_matched():
    pos = load_fasta_seqs(RAW / "q1_reviewed_positives.fasta")
    neg = load_fasta_seqs(RAW / "q1_matched_negatives.fasta")
    pairs, dropped = match_exact(pos, neg, widen=False)
    used_fallback = False
    if len(pairs) < MIN_PAIRS:
        pairs, dropped = match_exact(pos, neg, widen=True)
        used_fallback = True
    seqs = [pos[i] for i, _ in pairs] + [neg[j] for _, j in pairs]
    y = np.array([1]*len(pairs) + [0]*len(pairs))
    meta = {"pairs": len(pairs), "dropped_positives": dropped,
            "unused_negatives": len(neg) - len(pairs),
            "fallback_pm1_used": used_fallback}
    return seqs, y, meta

def make_splits(n):
    rng = np.random.RandomState(11)
    perm = rng.permutation(n)
    return [(np.concatenate([perm[:i*n//FOLDS], perm[(i+1)*n//FOLDS:]]),
             perm[i*n//FOLDS:(i+1)*n//FOLDS]) for i in range(FOLDS)]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", type=str, default="1-10")
    ap.add_argument("--falsifier", action="store_true")
    ap.add_argument("--finalize", action="store_true")
    args = ap.parse_args()
    seqs, y, meta = load_matched()
    if args.finalize:
        rows = [json.load(open(p)) for p in sorted(PART.glob("fold_*.json"),
                key=lambda p: int(p.name.split('_')[1].split('.')[0]))]
        frows = [json.load(open(p)) for p in sorted(FPART.glob("fold_*.json"),
                key=lambda p: int(p.name.split('_')[1].split('.')[0]))]
        assert len(rows) == FOLDS and len(frows) == FOLDS
        out = {"study": "study28_q7_length_confound", "matching": meta,
               "protocol": "study19 ensemble, rng(11) random 10-fold, exact-length 1:1 subset",
               "comparator_study22": {"auc": 0.9633, "mcc": 0.803},
               "folds": rows}
        for k in ["acc", "sens", "spec", "mcc", "auc"]:
            v = [r[k] for r in rows]; m, s = st.mean(v), st.stdev(v)
            out[k] = {"folds": [round(x,4) for x in v], "mean": round(m,4), "sd": round(s,4),
                      "lower95_one_sided": round(m - 1.833*s/math.sqrt(FOLDS), 4)}
        fmcc = st.mean(r["mcc"] for r in frows)
        out["falsifier_shuffle_mean_mcc"] = round(fmcc, 4)
        out["falsifier_valid"] = bool(abs(fmcc) < 0.1)
        out["delta_vs_study22_auc"] = round(out["auc"]["mean"] - 0.9633, 4)
        out["delta_vs_study22_mcc"] = round(out["mcc"]["mean"] - 0.803, 4)
        json.dump(out, open(FINAL, "w"), indent=1)
        print("FINAL:", json.dumps({"pairs": meta["pairs"],
              "mcc": out["mcc"]["mean"], "auc": out["auc"]["mean"],
              "delta_mcc": out["delta_vs_study22_mcc"], "falsifier_mcc": round(fmcc,4)}))
        return
    lo, hi = (int(x) for x in args.folds.split("-"))
    part = FPART if args.falsifier else PART
    part.mkdir(exist_ok=True)
    yy = y.copy()
    if args.falsifier:
        np.random.RandomState(11).shuffle(yy)
    splits = make_splits(len(seqs))
    for f in range(lo, hi+1):
        pf = part / f"fold_{f}.json"
        if pf.exists():
            print(f"fold {f}: exists, skip", flush=True); continue
        row = run_fold(f, seqs, yy, splits)
        json.dump(row, open(pf, "w"))
        print(f"fold {f}: acc {row['acc']:.4f} mcc {row['mcc']:.4f} auc {row['auc']:.4f} ({row['seconds']:.1f}s)", flush=True)

if __name__ == "__main__":
    main()
