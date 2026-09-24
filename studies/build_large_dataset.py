"""Build the scale-up AMP dataset with homology-aware splits.

Union-find over shared exact 10-mers links putative homologs in O(n * len),
then whole clusters are assigned to train/val/test so no near-duplicate
crosses splits. Negatives: reviewed non-AMP pool, length-bucket matched 1:1.
"""
import json, os, random, csv
from collections import defaultdict

ROOT = os.path.join(os.path.dirname(__file__), "..")
RAW = os.path.join(ROOT, "data", "raw")
PROC = os.path.join(ROOT, "data", "processed_large")
os.makedirs(PROC, exist_ok=True)
AA = set("ACDEFGHIKLMNPQRSTVWY")

def parse_fasta(path):
    recs, header, seq = [], None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if header is not None: recs.append((header, "".join(seq)))
            header, seq = line[1:], []
        elif line: seq.append(line)
    if header is not None: recs.append((header, "".join(seq)))
    return recs

def clean(recs):
    out, seen = [], set()
    for h, s in recs:
        s = s.upper()
        if not (10 <= len(s) <= 150) or not set(s) <= AA or s in seen: continue
        seen.add(s); out.append((h, s))
    return out

class DSU:
    def __init__(self, n): self.p = list(range(n))
    def find(self, x):
        while self.p[x] != x: self.p[x] = self.p[self.p[x]]; x = self.p[x]
        return x
    def union(self, a, b): self.p[self.find(a)] = self.find(b)

def main(seed=27, kmer=10):
    random.seed(seed)
    pos = clean(parse_fasta(os.path.join(RAW, "uniprot_kw0929_all_10-150.fasta")))
    neg_pool = clean(parse_fasta(os.path.join(RAW, "uniprot_nonamp_pool_10-150.fasta")))
    print(f"clean positives: {len(pos)}, clean negative pool: {len(neg_pool)}", flush=True)
    random.shuffle(neg_pool)
    pos_hist = defaultdict(int)
    for _, s in pos: pos_hist[len(s)//10] += 1
    neg, bc = [], defaultdict(int)
    for h, s in neg_pool:
        b = len(s)//10
        if bc[b] < pos_hist[b]:
            neg.append((h, s)); bc[b] += 1
        if len(neg) >= len(pos): break
    print(f"length-matched negatives: {len(neg)}", flush=True)
    allrec = [(1, h, s) for h, s in pos] + [(0, h, s) for h, s in neg]
    n = len(allrec)
    dsu = DSU(n)
    seen_kmer = {}
    for i, (_, _, s) in enumerate(allrec):
        for j in range(len(s) - kmer + 1):
            k = s[j:j+kmer]
            if k in seen_kmer: dsu.union(i, seen_kmer[k])
            else: seen_kmer[k] = i
    clusters = defaultdict(list)
    for i in range(n): clusters[dsu.find(i)].append(i)
    cids = list(clusters.keys()); random.shuffle(cids)
    nc = len(cids)
    tr = set(cids[:int(0.8*nc)]); va = set(cids[int(0.8*nc):int(0.9*nc)])
    counts = {"train": [0,0], "val": [0,0], "test": [0,0]}
    files = {sp: open(os.path.join(PROC, f"{sp}.csv"), "w", newline="") for sp in counts}
    writers = {}
    for sp, f in files.items():
        w = csv.writer(f); w.writerow(["label","header","sequence","length","cluster"]); writers[sp] = w
    for cid in cids:
        sp = "train" if cid in tr else "val" if cid in va else "test"
        for i in clusters[cid]:
            label, h, s = allrec[i]
            writers[sp].writerow([label, h, s, len(s), cid])
            counts[sp][0] += 1; counts[sp][1] += label
    for f in files.values(): f.close()
    stats = {"kmer": kmer, "seed": seed, "n_clean_pos": len(pos), "n_neg_used": len(neg),
             "n_records": n, "n_clusters": nc,
             "splits": {sp: {"n": c[0], "pos": c[1]} for sp, (c) in counts.items()}}
    json.dump(stats, open(os.path.join(PROC, "stats.json"), "w"), indent=2)
    print(json.dumps(stats, indent=2), flush=True)

if __name__ == "__main__":
    main()
