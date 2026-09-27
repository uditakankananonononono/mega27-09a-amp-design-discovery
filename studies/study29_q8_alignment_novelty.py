"""Study 29 - Q8 alignment-based novelty metric (addendum 8, locked 2026-09-27
BEFORE scoring). MMseqs2 easy-search (-s 7.5, --max-seqs 10000, e <= 1e-3, no
min-seq-id filter) of the 10 named InstiAMP-21 + 9 prior study12 candidates
(post-rerun 132df1a) vs the three pinned A2.2 corpora (APD6 3,306 / DRAMP
12,816 / DBAASP corrected catalog 16,239 sha 735a966c). Locked novelty flag:
NOVEL if no hit with pident >= 30% AND qcov >= 0.5 at e <= 1e-3 in any corpus.
Reported with the 3-mer Jaccard tiers as-is; neither lens overrides the other.
"""
import json, pathlib, subprocess, tempfile

ROOT = pathlib.Path(__file__).resolve().parent.parent
MMSEQS = "/tmp/mmseqs/bin/mmseqs"
OUT = ROOT / "results/study29_q8_alignment_novelty.json"
AA20 = set("ACDEFGHIKLMNPQRSTVWY")

def load_fasta(p):
    seqs, cur = [], None
    for line in open(p):
        line = line.strip()
        if line.startswith(">"):
            cur = ""
        elif line and cur is not None:
            cur += line
            seqs.append(cur)
    return seqs

def main():
    rescreen = json.load(open(ROOT / "results/study21_rescreen.json"))
    queries = [{"name": c["name"], "sequence": c["sequence"],
                "jaccard_tier": c["novelty"]["tier"]}
               for c in rescreen["named_candidates"] + rescreen["prior_candidates_rescored"]]
    dr = []
    for i, line in enumerate(open(ROOT / "data/raw/novelty_refs/dramp_natural_amps.txt")):
        if i == 0: continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) > 1 and parts[1] and set(parts[1]) <= AA20 and 5 <= len(parts[1]) <= 150:
            dr.append(parts[1])
    corpora = {"APD6": load_fasta(ROOT / "data/raw/apd6_natural_2024.fasta"),
               "DRAMP": dr,
               "DBAASP": load_fasta(ROOT / "data/raw/novelty_refs/dbaasp_all.fasta")}
    print({k: len(v) for k, v in corpora.items()}, flush=True)
    with tempfile.TemporaryDirectory() as td:
        qf = pathlib.Path(td) / "q.fa"
        qf.write_text("".join(f">{i}\n{q['sequence']}\n" for i, q in enumerate(queries)))
        results = {}
        for name, corpus in corpora.items():
            tf = pathlib.Path(td) / f"t_{name}.fa"
            tf.write_text("".join(f">{j}\n{s}\n" for j, s in enumerate(corpus)))
            m8 = pathlib.Path(td) / f"res_{name}.m8"
            subprocess.run([MMSEQS, "easy-search", str(qf), str(tf), str(m8), str(pathlib.Path(td)/f"tmp_{name}"),
                            "-s", "7.5", "--max-seqs", "10000", "-e", "1e-3", "--threads", "1"],
                           check=True, capture_output=True)
            best = {}
            for line in open(m8):
                f = line.rstrip("\n").split("\t")
                qi, ti, pid, ev = int(f[0]), int(f[1]), float(f[2]), float(f[10])
                qcov = (int(f[7]) - int(f[6]) + 1) / len(queries[qi]["sequence"])
                cur = best.get(qi)
                key = (pid >= 30 and qcov >= 0.5, pid, qcov, -ev)
                if cur is None or key > cur[0]:
                    best[qi] = (key, {"target_index": ti, "pident": round(pid, 1),
                                      "evalue": ev, "qcov": round(qcov, 3),
                                      "target_seq": corpus[ti]})
            results[name] = best
    out = {"study": "study29_q8_alignment_novelty",
           "engine": "MMseqs2 easy-search -s 7.5 --max-seqs 10000 -e 1e-3 (--threads 1, mechanics only)",
           "locked_flag": "NOVEL if no hit with pident>=30% AND qcov>=0.5 at e<=1e-3 in any corpus",
           "candidates": []}
    for i, q in enumerate(queries):
        per = {}
        novel = True
        for name in corpora:
            b = results[name].get(i)
            hit = b[1] if b else None
            flagged = bool(hit and hit["pident"] >= 30 and hit["qcov"] >= 0.5)
            if flagged: novel = False
            per[name] = {"best": hit, "flagged": flagged}
        out["candidates"].append({"name": q["name"], "sequence": q["sequence"],
                                  "jaccard_tier": q["jaccard_tier"],
                                  "per_corpus": per, "alignment_novel": novel})
    n_novel = sum(1 for c in out["candidates"] if c["alignment_novel"])
    out["summary"] = {"n": len(out["candidates"]), "alignment_novel": n_novel,
                      "jaccard_T1": sum(1 for c in out["candidates"] if c["jaccard_tier"] == "T1")}
    json.dump(out, open(OUT, "w"), indent=1)
    print(json.dumps(out["summary"]), flush=True)
    for c in out["candidates"]:
        print(c["name"], c["jaccard_tier"], "NOVEL" if c["alignment_novel"] else "HIT",
              {k: (v["best"]["pident"] if v["best"] else None) for k, v in c["per_corpus"].items()}, flush=True)

if __name__ == "__main__":
    main()
