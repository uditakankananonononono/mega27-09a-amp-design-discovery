"""Study 24 - Q4 per-family performance report (addendum-4 Q4 = A7).

No per-sequence prediction dumps exist in any prior study (all JSONs are
aggregate-only - disclosed). This reruns the EXACT study22 Q1 protocol
(same clean 2,973+2,973 reviewed dataset, same study19 ensemble, same
rng(11) random 10-fold splits, same threshold 0.5) but dumps per-sequence
out-of-fold ensemble probabilities, then reports metrics by motif-derived
family. Family rules (sequence-derived, locked in addendum-4 wording
"anionic, disulfide-rich, etc. by motif/annotation"):
  disulfide_rich : Cys fraction >= 0.10
  anionic        : net charge (D,E - K,R) < 0
  cationic       : net charge >= +2
  gly_rich       : Gly fraction >= 0.15
  other          : none of the above
Reported as-is per family: n, AUC (when both classes present), MCC, acc.
Sanity gate: pooled OOF metrics must reproduce study22 headline within
tolerance (seed-identical splits; torch nondeterminism may cause small
wobble - disclosed if observed).

  PYTHONPATH=. python3 studies/study24_q4_per_family.py --folds 1-10
  PYTHONPATH=. python3 studies/study24_q4_per_family.py --finalize
"""
from __future__ import annotations
import argparse, json, pathlib, time
import numpy as np
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score, accuracy_score, matthews_corrcoef

from pepdesign.models.cnn import PepCNN
from pepdesign.models.gnn import PepGNN
from pepdesign.models.ensemble import (train_binary, predict_scores,
                                       train_binary_shared_adj,
                                       predict_scores_shared_adj)
from studies.study14_feb2020_full import rich_features, load_fasta_seqs
from studies.study19_feb2020_v2 import (MAX_LEN, WEIGHTS, rich_features_v2,
                                        cnn_inputs)
from studies.study22_q1_reviewed_rebuild import load_data, make_splits

ROOT = pathlib.Path(__file__).resolve().parent.parent
PART = ROOT / "results/study24_q4_partial"
FINAL = ROOT / "results/study24_q4_per_family.json"
FOLDS = 10

def family(seq):
    n = len(seq)
    cys = seq.count("C") / n
    charge = seq.count("K") + seq.count("R") - seq.count("D") - seq.count("E")
    gly = seq.count("G") / n
    if cys >= 0.10:
        return "disulfide_rich"
    if charge < 0:
        return "anionic"
    if charge >= 2:
        return "cationic"
    if gly >= 0.15:
        return "gly_rich"
    return "other"

def run_fold_dump(fold_no, seqs, y, splits):
    tri, tei = splits[fold_no-1]
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
    return {"fold": fold_no, "tei": [int(i) for i in tei],
            "ens": [round(float(p), 5) for p in ens]}

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
        n = len(seqs)
        oof_p = np.zeros(n); seen = np.zeros(n, bool)
        for r in rows:
            for i, p in zip(r["tei"], r["ens"]):
                oof_p[i] = p; seen[i] = True
        assert seen.all()
        pred = (oof_p >= 0.5).astype(int)
        pooled = {"acc": float(accuracy_score(y, pred)),
                  "mcc": float(matthews_corrcoef(y, pred)),
                  "auc": float(roc_auc_score(y, oof_p))}
        fams = {}
        for i, s in enumerate(seqs):
            fams.setdefault(family(s), []).append(i)
        per = {}
        for fam, idx in sorted(fams.items()):
            fy = y[idx]; fp = oof_p[idx]; fb = pred[idx]
            rec = {"n": len(idx), "n_pos": int(fy.sum()),
                   "acc": float(accuracy_score(fy, fb)),
                   "mcc": float(matthews_corrcoef(fy, fb))}
            if fy.min() != fy.max():
                rec["auc"] = float(roc_auc_score(fy, fp))
            per[fam] = rec
        out = {"study": "study24_q4_per_family",
               "protocol": "study22 Q1 exact rerun with per-sequence OOF dumps",
               "pooled_oof": pooled,
               "study22_reference": {"acc": 0.9008, "mcc": 0.8030, "auc": 0.9633},
               "per_family": per}
        json.dump(out, open(FINAL, "w"), indent=1)
        print("FINAL:", json.dumps(out["pooled_oof"]), flush=True)
        print(json.dumps(per, indent=1), flush=True)
        return
    lo, hi = (int(x) for x in args.folds.split("-"))
    PART.mkdir(exist_ok=True)
    splits = make_splits(len(seqs))
    for f in range(lo, hi+1):
        pf = PART / f"fold_{f}.json"
        if pf.exists():
            print(f"fold {f}: exists, skip", flush=True); continue
        row = run_fold_dump(f, seqs, y, splits)
        json.dump(row, open(pf, "w"))
        print(f"fold {f}: dumped {len(row['tei'])} OOF predictions", flush=True)

if __name__ == "__main__":
    main()
