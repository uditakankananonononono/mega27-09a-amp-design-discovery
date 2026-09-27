"""Study 27 - Q3 DBAASP potency regression (addendum 6, locked 2026-09-27
before any label inspection; DEVIATION-1 + CORRECTION-1 disclosed in-doc).

Locked protocol (docs/PREREG_ADDENDUM6_2026-09-27.md):
- Assay filter: targetActivities with activityMeasureGroup.name == "MIC" only.
- Primary cohort: targetSpecies.name == "Escherichia coli ATCC 25922",
  unit ug/ml only; uM converted ONLY for plain 20-aa sequences
  (ug/ml = uM x MW/1000, average residue masses). Conversion count reported.
- Multimer records (chains under monomers[]) EXCLUDED - an MIC on a multimer
  is not an MIC on a chain (count disclosed).
- Replicates: median MIC per peptide. Target y = log10(median MIC ug/ml).
- Model: HistGradientBoostingRegressor(max_iter=300, random_state=11) on
  rich_features_v2 (study19 feature set).
- Splits: 10-fold cluster-held-out CV; MMseqs2 clusters (min-seq-id 0.3,
  cov 0.8, cov-mode 0, union-find single-linkage), rng(11) greedy assignment
  - identical machinery to study23.
- Null: global-median predictor on same folds.
- Falsifier: one label shuffle (rng 11); rho > 0.2 invalidates the arm.
- Metrics: Spearman rho (primary), MAE (log10 ug/ml), R^2, with null deltas.
  Reported as-is, no success threshold.
- Secondary (exploratory, labeled, no headline): all-species pooled ug/ml.

  PYTHONPATH=. python3 studies/study27_q3_potency.py --build
  PYTHONPATH=. python3 studies/study27_q3_potency.py --run
"""
import json, pathlib, sys
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
API = ROOT / "data/raw/dbaasp_api"
DS = ROOT / "results/study27_q3_dataset.json"
OUT = ROOT / "results/study27_q3_potency.json"

AA_MASS = {"A": 89.09, "R": 174.20, "N": 132.12, "D": 133.10, "C": 121.16,
           "E": 147.13, "Q": 146.15, "G": 75.07, "H": 155.16, "I": 131.17,
           "L": 131.17, "K": 146.19, "M": 149.21, "F": 165.19, "P": 115.13,
           "S": 105.09, "T": 119.12, "W": 204.23, "Y": 181.19, "V": 117.15}
AA20 = set(AA_MASS)
ECOLI = "Escherichia coli ATCC 25922"

def mw(seq):
    return sum(AA_MASS[a] for a in seq) + 18.02

def build():
    rows, n_multi, n_um = [], 0, 0
    recs = 0
    for p in sorted(API.glob("[0-9]*.json"), key=lambda x: int(x.stem)):
        d = json.loads(p.read_text())
        recs += 1
        seq = d.get("sequence")
        if not seq:
            n_multi += 1
            continue
        mics = []
        for ta in d.get("targetActivities", []):
            if (ta.get("activityMeasureGroup") or {}).get("name") != "MIC":
                continue
            sp = (ta.get("targetSpecies") or {}).get("name", "")
            unit = (ta.get("unit") or {}).get("name", "")
            try:
                val = float(ta.get("activity"))
            except (TypeError, ValueError):
                continue
            if unit == "µg/ml":
                ugml = val
            elif unit == "µM" and set(seq) <= AA20:
                ugml = val * mw(seq) / 1000.0
                n_um += 1
            else:
                continue
            mics.append((sp, ugml))
        if mics:
            rows.append({"id": d.get("id"), "sequence": seq, "assays": mics})
    def cohort(rows, species):
        out = []
        for r in rows:
            vals = [v for sp, v in r["assays"] if species is None or sp == species]
            if vals and 5 <= len(r["sequence"]) <= 150 and set(r["sequence"]) <= AA20:
                out.append({"id": r["id"], "sequence": r["sequence"],
                            "median_mic_ugml": float(np.median(vals)), "n_assays": len(vals)})
        return out
    primary = cohort(rows, ECOLI)
    secondary = cohort(rows, None)
    json.dump({"protocol": "addendum 6 locked build",
               "records_scanned": recs, "multimer_excluded": n_multi,
               "um_conversions": n_um,
               "primary_ecoli_atcc25922": primary,
               "secondary_all_species": secondary},
              open(DS, "w"))
    print(f"records={recs} multimers={n_multi} um_conv={n_um} "
          f"primary={len(primary)} secondary={len(secondary)}", flush=True)

def run():
    from scipy.stats import spearmanr
    from sklearn.ensemble import HistGradientBoostingRegressor
    from sklearn.metrics import mean_absolute_error, r2_score
    from studies.study19_feb2020_v2 import rich_features_v2
    from studies.study23_mmseqs2_split import mmseqs_clusters, make_splits
    ds = json.load(open(DS))
    pri = ds["primary_ecoli_atcc25922"]
    seqs = [r["sequence"] for r in pri]
    y = np.array([np.log10(r["median_mic_ugml"]) for r in pri])
    print(f"primary cohort n={len(seqs)}", flush=True)
    clusters = mmseqs_clusters(seqs)
    ncl = len(clusters)
    print(f"clusters={ncl}", flush=True)
    splits = make_splits(clusters, len(seqs))
    X = np.array([rich_features_v2(s) for s in seqs])
    def cv_pred(yy):
        pred = np.zeros(len(yy))
        for tri, tei in splits:
            m = HistGradientBoostingRegressor(max_iter=300, random_state=11)
            m.fit(X[tri], yy[tri])
            pred[tei] = m.predict(X[tei])
        return pred
    pred = cv_pred(y)
    null_pred = np.zeros(len(y))
    for tri, tei in splits:
        null_pred[tei] = np.median(y[tri])
    ysh = y.copy(); np.random.RandomState(11).shuffle(ysh)
    fals_pred = cv_pred(ysh)
    def metrics(yy, pp):
        return {"spearman": float(spearmanr(yy, pp).statistic),
                "mae": float(mean_absolute_error(yy, pp)),
                "r2": float(r2_score(yy, pp))}
    out = {"study": "study27_q3_potency", "cohort": "E. coli ATCC 25922, log10 median MIC ug/ml",
           "n": len(seqs), "n_clusters": ncl,
           "model": metrics(y, pred), "null_median": metrics(y, null_pred),
           "falsifier_shuffle": metrics(y, fals_pred),
           "falsifier_invalidates": bool(spearmanr(y, fals_pred).statistic > 0.2)}
    json.dump(out, open(OUT, "w"), indent=1)
    print(json.dumps(out, indent=1), flush=True)

if __name__ == "__main__":
    if "--build" in sys.argv:
        build()
    elif "--run" in sys.argv:
        run()
