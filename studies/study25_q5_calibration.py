"""Study 25 - Q5 calibration (addendum-4 Q5 = A8, locked 2026-09-27):
Platt + isotonic calibration on held-out folds; report Brier + slope.

Uses the study24 per-sequence out-of-fold ensemble probabilities (exact
study22 protocol rerun, 5,946 points). Protocol: 10-fold CV over the OOF
points (rng(11) assignment); for each fold fit calibrators on the other 9,
evaluate on the held-out fold; pool held-out predictions. Reported as-is:
uncalibrated vs Platt vs isotonic Brier score and calibration slope
(logistic recalibration slope of outcome on held-out calibrated scores).
No success threshold (per addendum-4 cheap-item wording).
"""
import json, pathlib
import numpy as np
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss

ROOT = pathlib.Path(__file__).resolve().parent.parent
PART = ROOT / "results/study24_q4_partial"
OUT = ROOT / "results/study25_q5_calibration.json"

from studies.study22_q1_reviewed_rebuild import load_data
seqs, y = load_data()
n = len(seqs)
oof = np.zeros(n)
for p in sorted(PART.glob("fold_*.json")):
    r = json.load(open(p))
    for i, pr in zip(r["tei"], r["ens"]):
        oof[i] = pr

rng = np.random.RandomState(11)
perm = rng.permutation(n)
folds = [perm[i*n//10:(i+1)*n//10] for i in range(10)]
eps = 1e-7
platt_p = np.zeros(n); iso_p = np.zeros(n)
for f in range(10):
    te = folds[f]
    tr = np.concatenate([folds[j] for j in range(10) if j != f])
    lr = LogisticRegression(max_iter=1000)
    lr.fit(oof[tr].reshape(-1, 1), y[tr])
    platt_p[te] = lr.predict_proba(oof[te].reshape(-1, 1))[:, 1]
    iso = IsotonicRegression(out_of_bounds="clip")
    iso.fit(oof[tr], y[tr])
    iso_p[te] = iso.predict(oof[te])

def slope(scores):
    lr = LogisticRegression(max_iter=1000)
    lr.fit(np.clip(scores, eps, 1-eps).reshape(-1, 1), y)
    return float(lr.coef_[0][0])

out = {"study": "study25_q5_calibration",
       "protocol": "10-fold held-out calibration over study24 OOF probabilities",
       "n": n,
       "uncalibrated": {"brier": float(brier_score_loss(y, oof)), "slope": slope(oof)},
       "platt": {"brier": float(brier_score_loss(y, platt_p)), "slope": slope(platt_p)},
       "isotonic": {"brier": float(brier_score_loss(y, iso_p)), "slope": slope(iso_p)}}
json.dump(out, open(OUT, "w"), indent=1)
print(json.dumps(out, indent=1))
