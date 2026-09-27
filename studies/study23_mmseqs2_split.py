"""Study 23 - Q2 MMseqs2 alignment-based clustering arm (addendum-4 Q2,
locked 2026-09-27): replaces the 3-mer Jaccard control with MMseqs2
alignment-based clusters (locked params: --min-seq-id 0.3 -c 0.8,
--cluster-mode 1 single-linkage/connected-component). Same dataset
(feb2020 2021+2021), same rng(11) greedy fold assignment, same 10-fold
cluster-held-out protocol, identical study19 ensemble. Comparator of
record: study21 t0.40 primary (mean MCC 0.8337). Reported as-is.

  python studies/study23_mmseqs2_split.py --clusters-only   # run MMseqs2
  python studies/study23_mmseqs2_split.py --folds 1-10
  python studies/study23_mmseqs2_split.py --finalize
"""
from __future__ import annotations
import argparse, json, pathlib, subprocess, tempfile

import numpy as np
import scipy.stats as sstats

from studies.study14_feb2020_full import load_fasta_seqs
from studies.study19_feb2020_v2 import run_fold

ROOT = pathlib.Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw/feb2020"
FOLDS = 10
PART = ROOT / "results/study23_partial"
MMSEQS = "/tmp/mmseqs/bin/mmseqs"
PARAMS = dict(min_seq_id=0.3, cov=0.8, cluster_mode=1)  # locked addendum-4 Q2

def mmseqs_clusters(seqs):
    with tempfile.TemporaryDirectory() as td:
        fa = pathlib.Path(td)/"all.fa"
        fa.write_text("".join(f">{i}\n{s}\n" for i, s in enumerate(seqs)))
        # align2clust in this MMseqs2 build rejects --cluster-mode 1, so
        # single-linkage is computed exactly: all-vs-all easy-search, keep
        # edges with pident >= 30% and bidirectional coverage >= 0.8
        # (cov-mode 0, the locked -c 0.8 semantics), then union-find.
        subprocess.run([MMSEQS, "easy-search", str(fa), str(fa), f"{td}/res.m8", f"{td}/tmp",
                        "--min-seq-id", str(PARAMS['min_seq_id']),
                        "-c", str(PARAMS['cov']), "--cov-mode", "0",
                        "-s", "7.5", "--max-seqs", "10000"],
                       check=True, capture_output=True)
        parent = list(range(len(seqs)))
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
        groups = {}
        for i in range(len(seqs)):
            groups.setdefault(find(i), []).append(i)
    return list(groups.values())

def make_splits(clusters, n):
    rng = np.random.default_rng(11)
    order = rng.permutation(len(clusters))
    folds = [[] for _ in range(FOLDS)]
    for ci in order:
        target = min(range(FOLDS), key=lambda f: len(folds[f]))
        folds[target].extend(clusters[ci])
    splits = []
    allset = set(range(n))
    for f in range(FOLDS):
        tei = np.array(sorted(folds[f]))
        tri = np.array(sorted(allset - set(folds[f])))
        splits.append((tri, tei))
    return splits

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", type=str, default=None)
    ap.add_argument("--finalize", action="store_true")
    ap.add_argument("--clusters-only", action="store_true")
    args = ap.parse_args()
    amp = load_fasta_seqs(RAW/"AMPS_02182020.fasta")
    dec = load_fasta_seqs(RAW/"DECOYS_02182020.fasta")
    seqs = amp + dec
    y = np.array([1]*len(amp) + [0]*len(dec))
    assert (len(amp), len(dec)) == (2021, 2021)
    cache = pathlib.Path("/tmp/study23_splits_cache.json")
    if cache.exists():
        d = json.loads(cache.read_text())
        clusters, splits = d["clusters"], [(np.array(a, dtype=int), np.array(b, dtype=int)) for a, b in d["splits"]]
        print("splits loaded from cache", flush=True)
    else:
        clusters = mmseqs_clusters(seqs)
        splits = make_splits(clusters, len(seqs))
        cache.write_text(json.dumps({"clusters": clusters,
            "splits": [(a.tolist(), b.tolist()) for a, b in splits]}))
    print(f"clusters: {len(clusters)} (largest {max(len(c) for c in clusters)}, singletons {sum(1 for c in clusters if len(c)==1)})", flush=True)
    (ROOT/"results/study23_clusters.json").write_text(json.dumps(
        {"tool": "MMseqs2 " + subprocess.run([MMSEQS, "version"], capture_output=True, text=True).stdout.strip(),
         "params": PARAMS, "n_clusters": len(clusters),
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
        out = {"protocol": "study23 locked (addendum4 Q2): cluster-held-out 10-fold CV, MMseqs2 min-seq-id 0.3 cov 0.8 single-linkage, rng(11) greedy assignment, identical study19 ensemble",
               "n_clusters": len(clusters),
               "study21_t040_primary_mean_mcc": 0.8337,
               "study19_random_split_mean_mcc": 0.8328,
               "summary": summary,
               "folds": rows}
        (ROOT/"results/study23_mmseqs2_split.json").write_text(json.dumps(out, indent=1))
        print(json.dumps({"summary": summary}, indent=1))
        return
    assert args.folds, "--folds N or N-M required"
    PART.mkdir(exist_ok=True)
    lo, hi = map(int, args.folds.split("-")) if "-" in args.folds else (int(args.folds), int(args.folds))
    for f in range(lo, hi+1):
        row = run_fold(f, seqs, y, splits)
        (PART/f"fold_{f}.json").write_text(json.dumps(row))
        print(f"fold {f}: acc {row['acc']:.4f} mcc {row['mcc']:.4f} auc {row['auc']:.4f}", flush=True)

if __name__ == "__main__":
    main()
