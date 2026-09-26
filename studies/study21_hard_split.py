"""Study 21 - HARDER family split (addendum 3, A3; locked 2026-09-27 BEFORE any
run). Identical protocol to study20 except 3-mer Jaccard threshold 0.40
(primary) with 0.35/0.45 sensitivity arms (--threshold). Locked criteria:
HARD-TEST SUCCESS mean MCC >= 0.799 (AMP Scanner published row); STRONG
SUCCESS mean MCC >= 0.8328 (study19 random-split mean); negative -> rule-6.

  python studies/study20_cluster_holdout.py --folds 1-10
  python studies/study20_cluster_holdout.py --finalize
"""
from __future__ import annotations
import argparse, json, pathlib, time
from collections import defaultdict

import numpy as np
import scipy.stats as sstats

from studies.study14_feb2020_full import load_fasta_seqs
from studies.study19_feb2020_v2 import run_fold

ROOT = pathlib.Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw/feb2020"
FOLDS = 10
JACCARD_THRESHOLD = float(__import__("sys").argv[__import__("sys").argv.index("--threshold")+1]) if "--threshold" in __import__("sys").argv else 0.40
PART = ROOT / f"results/study21_partial_t{JACCARD_THRESHOLD}"

def kmers(s, k=3):
    return {s[i:i+k] for i in range(len(s)-k+1)}

def cluster(seqs):
    sets = [kmers(s) for s in seqs]
    index = defaultdict(list)
    for i, st in enumerate(sets):
        for km in st:
            index[km].append(i)
    parent = list(range(len(seqs)))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb
    for i, st in enumerate(sets):
        cand = defaultdict(int)
        for km in st:
            for j in index[km]:
                if j > i:
                    cand[j] += 1
        for j, shared in cand.items():
            if shared >= JACCARD_THRESHOLD * (len(st) + len(sets[j])) / (1 + JACCARD_THRESHOLD):
                jac = shared / (len(st) + len(sets[j]) - shared)
                if jac >= JACCARD_THRESHOLD:
                    union(i, j)
        if i % 500 == 0:
            print(f"  cluster pass {i}/{len(seqs)}", flush=True)
    groups = defaultdict(list)
    for i in range(len(seqs)):
        groups[find(i)].append(i)
    return list(groups.values())

def make_splits(seqs, y):
    clusters = cluster(seqs)
    rng = np.random.default_rng(11)
    order = rng.permutation(len(clusters))
    folds = [[] for _ in range(FOLDS)]
    for ci in order:
        target = min(range(FOLDS), key=lambda f: len(folds[f]))
        folds[target].extend(clusters[ci])
    splits = []
    allset = set(range(len(seqs)))
    for f in range(FOLDS):
        tei = np.array(sorted(folds[f]))
        tri = np.array(sorted(allset - set(folds[f])))
        splits.append((tri, tei))
    return clusters, splits

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", type=str, default=None)
    ap.add_argument("--finalize", action="store_true")
    ap.add_argument("--clusters-only", action="store_true")
    ap.add_argument("--threshold", type=float, default=0.40)
    args = ap.parse_args()
    amp = load_fasta_seqs(RAW/"AMPS_02182020.fasta")
    dec = load_fasta_seqs(RAW/"DECOYS_02182020.fasta")
    seqs = amp + dec
    y = np.array([1]*len(amp) + [0]*len(dec))
    assert (len(amp), len(dec)) == (2021, 2021)
    cache = pathlib.Path("/tmp/study21_splits_cache_t%s.json" % JACCARD_THRESHOLD)
    if cache.exists():
        d = json.loads(cache.read_text())
        clusters, splits = d["clusters"], [(tuple(map(int, a)), tuple(map(int, b))) for a, b in d["splits"]]
        splits = [(np.array(tri), np.array(tei)) for tri, tei in splits]
        print("splits loaded from cache", flush=True)
    else:
        clusters, splits = make_splits(seqs, y)
        cache.write_text(json.dumps({"clusters": clusters,
            "splits": [(list(map(int, a)), list(map(int, b))) for a, b in splits]}))
    print(f"clusters: {len(clusters)} (largest {max(len(c) for c in clusters)}, singletons {sum(1 for c in clusters if len(c)==1)})", flush=True)
    (ROOT/f"results/study21_clusters_t{JACCARD_THRESHOLD}.json").write_text(json.dumps(
        {"n_clusters": len(clusters), "threshold_3mer_jaccard": JACCARD_THRESHOLD,
         "sizes": sorted((len(c) for c in clusters), reverse=True)[:50],
         "fold_sizes": [len(te) for _, te in splits]}))
    if args.clusters_only:
        return
    if args.finalize:
        rows = [json.loads(p.read_text()) for p in sorted(PART.glob("fold_*.json"))]
        assert len(rows) == FOLDS, f"only {len(rows)} folds present"
        summary = {}
        for m in ("sens","spec","acc","mcc","auc"):
            vals = np.array([r[m] for r in rows])
            t95 = sstats.t.ppf(0.95, len(vals)-1) * vals.std(ddof=1)/np.sqrt(len(vals))
            summary[m] = {"mean": float(vals.mean()), "sd": float(vals.std(ddof=1)),
                          "one_sided_95_lower": float(vals.mean()-t95)}
        out = {"protocol": "study21 locked (addendum3 A3): cluster-held-out 10-fold CV, single-linkage 3-mer Jaccard>=%s, rng(11) greedy assignment, identical study19 ensemble" % JACCARD_THRESHOLD,
               "n_clusters": len(clusters),
               "study19_random_split_mean_mcc": 0.8328,
               "summary": summary,
               "hard_test_success": bool(summary["mcc"]["mean"] >= 0.799),
               "strong_success": bool(summary["mcc"]["mean"] >= 0.8328),
               "folds": rows}
        (ROOT/f"results/study21_hard_split_t{JACCARD_THRESHOLD}.json").write_text(json.dumps(out, indent=1))
        print(json.dumps({"summary": summary, "hard_test_success": out["hard_test_success"], "strong_success": out["strong_success"]}, indent=1))
        return
    assert args.folds, "--folds N or N-M required"
    PART.mkdir(exist_ok=True)
    lo, hi = map(int, args.folds.split("-")) if "-" in args.folds else (int(args.folds), int(args.folds))
    for f in range(lo, hi+1):
        row = run_fold(f, seqs, y, splits)
        (PART/f"fold_{f}.json").write_text(json.dumps(row))
        print(f"fold {f}: acc {row['acc']:.4f} mcc {row['mcc']:.4f} auc {row['auc']:.4f} ({row['seconds']}s)", flush=True)

if __name__ == "__main__":
    main()
