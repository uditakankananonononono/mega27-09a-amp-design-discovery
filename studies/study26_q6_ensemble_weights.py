"""Study 26 - Q6 ensemble-weight arm (addendum-4 Q6 = A9) + branch-level OOF dumps.

Q6: ensemble weights optimized on validation folds vs simple average; locked
comparison, paired folds. Needs branch-level out-of-fold scores (study24
dumped ensemble only), so this reruns the identical study22 protocol once
more and dumps every branch score per test sequence. Weight optimization
and the locked comparison happen at finalize: on each fold pair, weights
optimized on the training-side OOF scores (Nelder-Mead on MCC, simplex
seed rng(11)) are compared against simple average and the locked study19
weights on the paired held-out fold. Reported as-is.

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
PART = ROOT / "results/study26_q6_partial"
FINAL = ROOT / "results/study26_q6_ensemble_weights.json"
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
            "ens": [round(float(p), 5) for p in ens],
            "branch": {k: [round(float(p), 5) for p in v] for k, v in branch.items()}}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", type=str, default="1-10")
    ap.add_argument("--finalize", action="store_true")
    args = ap.parse_args()
    seqs, y = load_data()
    if args.finalize:
        from scipy.optimize import minimize
        rows = [json.load(open(p)) for p in sorted(PART.glob("fold_*.json"),
                key=lambda p: int(p.name.split('_')[1].split('.')[0]))]
        assert len(rows) == FOLDS, len(rows)
        n = len(seqs)
        B = {k: np.zeros(n) for k in WEIGHTS}
        seen = np.zeros(n, bool)
        for r in rows:
            for j, i in enumerate(r["tei"]):
                for k in WEIGHTS:
                    B[k][i] = r["branch"][k][j]
                seen[i] = True
        assert seen.all()
        names = list(WEIGHTS)
        M = np.stack([B[k] for k in names])
        rng = np.random.RandomState(11)
        perm = rng.permutation(n)
        pf = [perm[i*n//FOLDS:(i+1)*n//FOLDS] for i in range(FOLDS)]
        def mcc_of(w, idx):
            s = w @ M[:, idx]
            return matthews_corrcoef(y[idx], (s >= 0.5).astype(int))
        res = {"locked_study19": [], "simple_average": [], "optimized": []}
        for f in range(FOLDS):
            te = pf[f]
            tr = np.concatenate([pf[j] for j in range(FOLDS) if j != f])
            w0 = np.array([WEIGHTS[k] for k in names])
            opt = minimize(lambda w: -mcc_of(w, tr), w0, method="Nelder-Mead",
                           options={"maxiter": 2000, "xatol": 1e-4, "fatol": 1e-6})
            res["locked_study19"].append(mcc_of(w0, te))
            res["simple_average"].append(mcc_of(np.ones(len(names))/len(names), te))
            res["optimized"].append(mcc_of(opt.x, te))
        import statistics as st
        out = {"study": "study26_q6_ensemble_weights",
               "protocol": "addendum-4 Q6: weights optimized on validation folds vs simple average; paired folds; reported as-is",
               "branches": names,
               "folds": res,
               "mean_mcc": {k: round(st.mean(v), 4) for k, v in res.items()},
               "deltas_vs_locked": {k: round(st.mean(v) - st.mean(res["locked_study19"]), 4)
                                    for k, v in res.items() if k != "locked_study19"}}
        json.dump(out, open(FINAL, "w"), indent=1)
        print("FINAL:", json.dumps(out["mean_mcc"]), json.dumps(out["deltas_vs_locked"]), flush=True)
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
        print(f"fold {f}: dumped {len(row['tei'])} OOF predictions (+branches)", flush=True)

if __name__ == "__main__":
    main()
