"""Study 19 - locked improvement on the Feb2020 benchmark (see
docs/PREREG_ADDENDUM_2026-09-26.md). Identical folds/base learners to study 14,
plus HGB branch on rich_features_v2 and locked a-priori ensemble weights
(CNN7 .20, CNN27 .20, GNN .15, RF .20, HGB .25), threshold 0.5.

Usage:
  python studies/study19_feb2020_v2.py --folds 1        # one fold -> partial json
  python studies/study19_feb2020_v2.py --finalize       # aggregate -> results/study19_feb2020_v2.json
"""
from __future__ import annotations
import argparse, json, pathlib, time

import numpy as np
import torch
from scipy import stats as sstats
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, accuracy_score, matthews_corrcoef

from pepdesign.encoding import combined_features
from pepdesign.models.cnn import PepCNN
from pepdesign.models.gnn import PepGNN
from pepdesign.models.ensemble import (train_binary, predict_scores,
                                       train_binary_shared_adj,
                                       predict_scores_shared_adj)
from studies.study14_feb2020_full import (rich_features, load_fasta_seqs,
                                          cnn_inputs, KD, BOMAN, AA)

ROOT = pathlib.Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw/feb2020"
PART = ROOT / "results/study19_partial"
MAX_LEN = 200
FOLDS = 10
CNN_SEEDS = (7, 27)
WEIGHTS = {"cnn7": 0.20, "cnn27": 0.20, "gnn": 0.15, "rf": 0.20, "hgb": 0.25}

ARO = set("FWY")
ALI = set("ILV")
PKA_C = {"C": 8.3, "Y": 10.1, "D": 3.9, "E": 4.1, "H": 6.0, "K": 10.5, "R": 12.4}
GROUPS = [set("AVLIM"), set("FYW"), set("STNQ"), set("KRH"), set("DE"), set("C"), set("G"), set("P")]

def est_pi(seq):
    def charge(ph):
        pos = 10**PKA_C["K"]/(10**ph+10**PKA_C["K"]) * seq.count("K") \
            + 10**PKA_C["R"]/(10**ph+10**PKA_C["R"]) * seq.count("R") \
            + 10**PKA_C["H"]/(10**ph+10**PKA_C["H"]) * seq.count("H") \
            + 10**7.0/(10**ph+10**7.0)
        neg = 10**ph/(10**ph+10**PKA_C["D"]) * seq.count("D") \
            + 10**ph/(10**ph+10**PKA_C["E"]) * seq.count("E") \
            + 10**ph/(10**ph+10**PKA_C["C"]) * seq.count("C") \
            + 10**ph/(10**ph+10**PKA_C["Y"]) * seq.count("Y") \
            + 10**ph/(10**ph+10**2.34)
        return pos - neg
    lo, hi = 0.0, 14.0
    for _ in range(40):
        mid = (lo+hi)/2
        if charge(mid) > 0: lo = mid
        else: hi = mid
    return (lo+hi)/2

def hydrophobic_moment(seq, angle=100.0, window=11):
    if len(seq) < window: window = len(seq)
    rad = np.deg2rad(angle); best = 0.0
    for i in range(len(seq)-window+1):
        w = seq[i:i+window]
        x = sum(KD.get(a,0.0)*np.cos(j*rad) for j,a in enumerate(w))
        y = sum(KD.get(a,0.0)*np.sin(j*rad) for j,a in enumerate(w))
        best = max(best, (x*x+y*y)**0.5/window)
    return best

def max_window_kd(seq, window=9):
    if len(seq) <= window:
        return float(np.mean([KD.get(a,0.0) for a in seq]))
    return max(float(np.mean([KD.get(a,0.0) for a in seq[i:i+window]]))
               for i in range(len(seq)-window+1))

def rich_features_v2(seq):
    v = list(rich_features(seq))
    n = len(seq)
    v.append(sum(1 for a in seq if a in ARO)/n)          # aromaticity
    v.append(sum(1 for a in seq if a in ALI)/n)          # aliphatic fraction
    v.append(est_pi(seq))                                # estimated pI
    v.append(max_window_kd(seq))                         # max hydrophobic stretch
    v.append(hydrophobic_moment(seq, angle=160.0))       # beta moment
    v += [sum(1 for a in seq if a in g)/n for g in GROUPS]  # grouped composition
    return v

def run_fold(fold_no, seqs, y, splits):
    tri, tei = splits[fold_no-1]
    t0 = time.time()
    tr_s = [seqs[i] for i in tri]; tr_y = y[tri]
    te_s = [seqs[i] for i in tei]; te_y = y[tei]
    branch = {}
    for name, seed in (("cnn7", 7), ("cnn27", 27)):
        cnn = train_binary(PepCNN(), cnn_inputs, tr_s, tr_y, MAX_LEN,
                           epochs=20, batch_size=64, pos_weight=1.0, seed=seed)
        branch[name] = predict_scores(cnn, cnn_inputs, te_s, MAX_LEN)
    gnn = train_binary_shared_adj(PepGNN(), tr_s, tr_y, MAX_LEN, epochs=10,
                                  batch_size=64, pos_weight=1.0, seed=11)
    branch["gnn"] = predict_scores_shared_adj(gnn, te_s, MAX_LEN)
    Xtr = np.array([rich_features(s) for s in tr_s])
    Xte = np.array([rich_features(s) for s in te_s])
    rf = RandomForestClassifier(n_estimators=300, n_jobs=2, random_state=42)
    rf.fit(Xtr, tr_y)
    branch["rf"] = rf.predict_proba(Xte)[:, 1]
    Xtr2 = np.array([rich_features_v2(s) for s in tr_s])
    Xte2 = np.array([rich_features_v2(s) for s in te_s])
    hgb = HistGradientBoostingClassifier(random_state=42)
    hgb.fit(Xtr2, tr_y)
    branch["hgb"] = hgb.predict_proba(Xte2)[:, 1]
    ens = sum(WEIGHTS[k]*branch[k] for k in WEIGHTS)
    pred = (ens >= 0.5).astype(int)
    row = {"fold": fold_no,
           "sens": float(((pred==1)&(te_y==1)).sum()/max(((te_y==1)).sum(),1)),
           "spec": float(((pred==0)&(te_y==0)).sum()/max(((te_y==0)).sum(),1)),
           "acc": float(accuracy_score(te_y, pred)),
           "mcc": float(matthews_corrcoef(te_y, pred)),
           "auc": float(roc_auc_score(te_y, ens)),
           "branch_auc": {k: float(roc_auc_score(te_y, v)) for k, v in branch.items()},
           "seconds": round(time.time()-t0, 1)}
    return row

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", type=str, default=None)
    ap.add_argument("--finalize", action="store_true")
    args = ap.parse_args()
    amp = load_fasta_seqs(RAW/"AMPS_02182020.fasta")
    dec = load_fasta_seqs(RAW/"DECOYS_02182020.fasta")
    seqs = amp + dec
    y = np.array([1]*len(amp) + [0]*len(dec))
    assert (len(amp), len(dec)) == (2021, 2021)
    skf = StratifiedKFold(n_splits=FOLDS, shuffle=True, random_state=42)
    splits = list(skf.split(seqs, y))
    if args.finalize:
        rows = [json.loads(p.read_text()) for p in sorted(PART.glob("fold_*.json"))]
        assert len(rows) == FOLDS, f"only {len(rows)} folds present"
        out = {"benchmark": "AMP Scanner Vr.2 Feb2020 dataset (dveltri.com)",
               "n_amp": len(amp), "n_decoy": len(dec),
               "protocol": "study19 locked: study14 folds + HGB branch on rich_features_v2, a-priori weights CNN7 .20/CNN27 .20/GNN .15/RF .20/HGB .25, thr 0.5",
               "weights": WEIGHTS,
               "published_feb2020_cv": {"sens":0.906,"spec":0.891,"acc":0.899,"mcc":0.799,"auc":0.962},
               "folds": rows}
        summary = {}
        for m in ("sens","spec","acc","mcc","auc"):
            vals = np.array([r[m] for r in rows])
            t95 = sstats.t.ppf(0.95, len(vals)-1) * vals.std(ddof=1)/np.sqrt(len(vals))
            summary[m] = {"mean": float(vals.mean()), "sd": float(vals.std(ddof=1)),
                          "one_sided_95_lower": float(vals.mean()-t95)}
        out["summary"] = summary
        pub = out["published_feb2020_cv"]
        out["gates"] = {
            "G1_mcc_beat": bool(summary["mcc"]["mean"] > pub["mcc"] and summary["mcc"]["one_sided_95_lower"] > pub["mcc"]),
            "G2_acc_beat": bool(summary["acc"]["mean"] > pub["acc"]),
            "G3_auc_beat": bool(summary["auc"]["mean"] > pub["auc"])}
        (ROOT/"results/study19_feb2020_v2.json").write_text(json.dumps(out, indent=1))
        print(json.dumps({"summary": summary, "gates": out["gates"]}, indent=1))
        return
    assert args.folds, "--folds N or N-M required"
    PART.mkdir(exist_ok=True)
    if "-" in args.folds:
        lo, hi = map(int, args.folds.split("-"))
        fold_list = list(range(lo, hi+1))
    else:
        fold_list = [int(args.folds)]
    for f in fold_list:
        row = run_fold(f, seqs, y, splits)
        (PART/f"fold_{f}.json").write_text(json.dumps(row))
        print(f"fold {f}: acc {row['acc']:.4f} mcc {row['mcc']:.4f} auc {row['auc']:.4f} ({row['seconds']}s)", flush=True)

if __name__ == "__main__":
    main()
