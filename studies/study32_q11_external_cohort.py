"""Study 32 - Q11 external cohort (addendum 9, locked 2026-09-27 BEFORE
scoring). DBAASP activity-defined E. coli ATCC 25922 cohort from the completed
blind pull: positives median MIC <= 32 ug/ml, negatives ALL reported MICs
> 256 ug/ml (uM converted per study27). Contamination control: exact dedup vs
all training seqs + MMseqs2 cluster-level removal (study23 params, cohort and
training clustered jointly; cohort seqs in training-touching clusters
dropped). Model: study22-recipe ensemble trained on the FULL Q1 dataset
(2,973+2,973, study21 final-model pattern, no threshold/weight changes).
Metrics: AUROC primary, MCC @ locked 0.5 secondary, distributions; Spearman vs
log10 MIC on positives as study27 context. Reported as-is, any direction.

  PYTHONPATH=. python3 studies/study32_q11_external_cohort.py --build
  PYTHONPATH=. python3 studies/study32_q11_external_cohort.py --run
"""
from __future__ import annotations
import json, pathlib, subprocess, sys, tempfile
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
API = ROOT / "data/raw/dbaasp_api"
COHORT = ROOT / "results/study32_q11_cohort.json"
OUT = ROOT / "results/study32_q11_external_cohort.json"
MMSEQS = "/tmp/mmseqs/bin/mmseqs"
ECOLI = "Escherichia coli ATCC 25922"

sys.path.insert(0, str(ROOT))
from studies.study27_q3_potency import AA20, mw  # locked conversion reuse
from studies.study14_feb2020_full import load_fasta_seqs

def build():
    rows = []
    for p in sorted(API.glob("[0-9]*.json"), key=lambda x: int(x.stem)):
        d = json.loads(p.read_text())
        seq = d.get("sequence")
        if not seq or not (5 <= len(seq) <= 150) or not set(seq) <= AA20:
            continue
        vals = []
        for ta in d.get("targetActivities", []):
            if (ta.get("activityMeasureGroup") or {}).get("name") != "MIC":
                continue
            if (ta.get("targetSpecies") or {}).get("name", "") != ECOLI:
                continue
            unit = (ta.get("unit") or {}).get("name", "")
            try:
                val = float(ta.get("activity"))
            except (TypeError, ValueError):
                continue
            if unit == "µg/ml":
                vals.append(val)
            elif unit == "µM":
                vals.append(val * mw(seq) / 1000.0)
        if vals:
            rows.append({"id": d.get("id"), "sequence": seq, "mics": vals})
    pos = [r for r in rows if float(np.median(r["mics"])) <= 32.0]
    neg = [r for r in rows if min(r["mics"]) > 256.0]
    train = set(load_fasta_seqs(ROOT/"data/raw/q1_reviewed_positives.fasta"))
    train |= set(load_fasta_seqs(ROOT/"data/raw/q1_matched_negatives.fasta"))
    train |= set(load_fasta_seqs(ROOT/"data/raw/feb2020/AMPS_02182020.fasta"))
    train |= set(load_fasta_seqs(ROOT/"data/raw/feb2020/DECOYS_02182020.fasta"))
    n_exact_pos = sum(1 for r in pos if r["sequence"] in train)
    n_exact_neg = sum(1 for r in neg if r["sequence"] in train)
    pos = [r for r in pos if r["sequence"] not in train]
    neg = [r for r in neg if r["sequence"] not in train]
    # joint clustering: training seqs + cohort seqs, study23 params
    tseqs = sorted(train)
    allseqs = tseqs + [r["sequence"] for r in pos] + [r["sequence"] for r in neg]
    with tempfile.TemporaryDirectory() as td:
        fa = pathlib.Path(td)/"all.fa"
        fa.write_text("".join(f">{i}\n{s}\n" for i, s in enumerate(allseqs)))
        subprocess.run([MMSEQS, "easy-search", str(fa), str(fa), f"{td}/res.m8", f"{td}/tmp",
                        "--min-seq-id", "0.3", "-c", "0.8", "--cov-mode", "0",
                        "-s", "7.5", "--max-seqs", "10000", "--threads", "1"],
                       check=True, capture_output=True)
        parent = list(range(len(allseqs)))
        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]; x = parent[x]
            return x
        for line in open(f"{td}/res.m8"):
            f = line.split("\t")
            q, m = int(f[0]), int(f[1])
            if q != m:
                rq, rm = find(q), find(m)
                if rq != rm: parent[rq] = rm
    tainted = set()
    groups = {}
    for i in range(len(allseqs)):
        groups.setdefault(find(i), []).append(i)
    for g in groups.values():
        if any(i < len(tseqs) for i in g):
            tainted.update(i for i in g if i >= len(tseqs))
    np_, nn_ = len(pos), len(neg)
    pos = [r for k, r in enumerate(pos) if (len(tseqs)+k) not in tainted]
    neg = [r for k, r in enumerate(neg) if (len(tseqs)+np_+k) not in tainted]
    meta = {"records_with_ecoli_mic": len(rows),
            "pos_pre": np_, "neg_pre": nn_,
            "exact_dup_dropped_pos": n_exact_pos, "exact_dup_dropped_neg": n_exact_neg,
            "pos_final": len(pos), "neg_final": len(neg)}
    meta["cluster_dropped_pos"] = np_ - n_exact_pos - len(pos)
    meta["cluster_dropped_neg"] = nn_ - n_exact_neg - len(neg)
    json.dump({"protocol": "addendum 9 locked build", "meta": meta,
               "positives": [{"sequence": r["sequence"], "median_mic_ugml": float(np.median(r["mics"]))} for r in pos],
               "negatives": [{"sequence": r["sequence"], "min_mic_ugml": float(min(r["mics"]))} for r in neg]},
              open(COHORT, "w"))
    print(json.dumps(meta), flush=True)

def run():
    import torch
    from studies.study14_feb2020_full import (rich_features, train_binary,
        train_binary_shared_adj, MAX_LEN)
    from studies.study19_feb2020_v2 import rich_features_v2
    from pepdesign.models.cnn import PepCNN
    from pepdesign.models.gnn import PepGNN
    from pepdesign.models.ensemble import predict_scores, predict_scores_shared_adj
    from pepdesign.encoding import combined_features
    from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier
    from sklearn.metrics import roc_auc_score, matthews_corrcoef
    from scipy.stats import spearmanr

    pos_tr = load_fasta_seqs(ROOT/"data/raw/q1_reviewed_positives.fasta")
    neg_tr = load_fasta_seqs(ROOT/"data/raw/q1_matched_negatives.fasta")
    seqs = pos_tr + neg_tr
    y = np.array([1]*len(pos_tr) + [0]*len(neg_tr))
    assert len(pos_tr) == len(neg_tr) == 2973
    print("training study22-recipe ensemble on full Q1 dataset ...", flush=True)
    def inputs(ss, ml):
        return (torch.stack([torch.from_numpy(combined_features(x, ml)) for x in ss]),)
    branches = {}
    for name, seed in (("cnn7", 7), ("cnn27", 27)):
        branches[name] = train_binary(PepCNN(), inputs, seqs, y, MAX_LEN,
                                      epochs=20, batch_size=64, pos_weight=1.0, seed=seed)
        print(f"  {name} trained", flush=True)
    branches["gnn"] = train_binary_shared_adj(PepGNN(), seqs, y, MAX_LEN, epochs=10,
                                              batch_size=64, pos_weight=1.0, seed=11)
    print("  gnn trained", flush=True)
    X = np.array([rich_features(s) for s in seqs])
    rf = RandomForestClassifier(n_estimators=300, n_jobs=2, random_state=42).fit(X, y)
    X2 = np.array([rich_features_v2(s) for s in seqs])
    hgb = HistGradientBoostingClassifier(random_state=42).fit(X2, y)
    print("  rf+hgb trained", flush=True)
    def ens(s):
        c7 = predict_scores(branches["cnn7"], inputs, [s], MAX_LEN)[0]
        c27 = predict_scores(branches["cnn27"], inputs, [s], MAX_LEN)[0]
        g = predict_scores_shared_adj(branches["gnn"], [s], MAX_LEN)[0]
        r = rf.predict_proba(np.array([rich_features(s)]))[:, 1][0]
        h = hgb.predict_proba(np.array([rich_features_v2(s)]))[:, 1][0]
        return float(.20*c7 + .20*c27 + .15*g + .20*r + .25*h)

    co = json.load(open(COHORT))
    pseq = [r["sequence"] for r in co["positives"]]
    nseq = [r["sequence"] for r in co["negatives"]]
    print(f"scoring external cohort: {len(pseq)} pos, {len(nseq)} neg ...", flush=True)
    ps = [ens(s) for s in pseq]
    ns = [ens(s) for s in nseq]
    yy = np.array([1]*len(ps) + [0]*len(ns))
    sc = np.array(ps + ns)
    pred = (sc >= 0.5).astype(int)
    out = {"study": "study32_q11_external_cohort", "cohort_meta": co["meta"],
           "model": "study22-recipe ensemble, full Q1 dataset, locked thr 0.5",
           "auroc": float(roc_auc_score(yy, sc)),
           "mcc_at_0.5": float(matthews_corrcoef(yy, pred)),
           "pos_scores": {"n": len(ps), "median": float(np.median(ps)),
                          "q1": float(np.percentile(ps, 25)), "q3": float(np.percentile(ps, 75))},
           "neg_scores": {"n": len(ns), "median": float(np.median(ns)),
                          "q1": float(np.percentile(ns, 25)), "q3": float(np.percentile(ns, 75))},
           "pos_spearman_vs_log_mic": float(spearmanr(
               [np.log10(r["median_mic_ugml"]) for r in co["positives"]], ps).statistic)}
    json.dump(out, open(OUT, "w"), indent=1)
    print(json.dumps(out, indent=1), flush=True)

if __name__ == "__main__":
    if "--build" in sys.argv:
        build()
    elif "--run" in sys.argv:
        run()
