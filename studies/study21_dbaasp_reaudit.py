"""study21 DBAASP-arm re-audit (addendum 6 DEVIATION-1 follow-up).

The original study21 novelty audit read the DEGENERATE dbaasp_all.fasta
(546 entries / 12 unique sequences). This re-audit recomputes ONLY the
DBAASP arm of every screened candidate against the corrected full-catalog
reference (dbaasp_catalog_2026-09-27.fasta), recomputes overall/tier, and
writes results/study21_dbaasp_reaudit.json. APD6/DRAMP arms are untouched.
"""
import json, pathlib
from studies.study12_candidates import kmer_jaccard, longest_common_substr
from studies.study14_feb2020_full import load_fasta_seqs

ROOT = pathlib.Path(__file__).resolve().parent.parent
CAT = ROOT / "data/raw/novelty_refs/dbaasp_catalog_2026-09-27.fasta"
SRC = ROOT / "results/study21_rescreen.json"
OUT = ROOT / "results/study21_dbaasp_reaudit.json"

def audit(s, corpus):
    jmax, lmax = 0.0, 0
    for r in corpus:
        j = kmer_jaccard(s, r)
        if j > jmax:
            jmax = j
        if jmax >= 0.35:
            l = longest_common_substr(s, r)
            if l > lmax:
                lmax = l
    return round(jmax, 3), lmax

def main():
    corpus = load_fasta_seqs(CAT)
    d = json.loads(SRC.read_text())
    groups = ["named_candidates", "prior_candidates_rescored", "rediscovery_controls"]
    changed, rows = [], {}
    for g in groups:
        out = []
        for c in d.get(g, []):
            s = c["sequence"]
            j, l = audit(s, corpus)
            old = [c["novelty"]["per_db"]["DBAASP"]["max_3mer_jaccard"],
                   c["novelty"]["per_db"]["DBAASP"]["lcs"]]
            rec = {"name": c.get("name") or c.get("id") or c.get("header", "")[:60],
                   "old_DBAASP": old, "new_DBAASP": [j, l]}
            # recompute overall/tier with the new DBAASP arm
            per = {k: [v["max_3mer_jaccard"], v["lcs"]] for k, v in c["novelty"]["per_db"].items()}
            per["DBAASP"] = [j, l]
            overall = max(v[0] for v in per.values())
            tier = "T1" if overall < 0.35 else ("T2" if overall <= 0.6 else "T3")
            rec["old_overall"] = c["novelty"].get("overall_max_jaccard")
            rec["new_overall"] = round(overall, 3)
            rec["old_tier"] = c["novelty"].get("tier")
            rec["new_tier"] = tier
            if rec["old_DBAASP"] != [j, l] or rec["old_tier"] != tier:
                changed.append(rec)
            out.append(rec)
        rows[g] = out
    OUT.write_text(json.dumps({
        "protocol": "DBAASP-arm re-audit vs corrected full-catalog reference "
                    "(addendum 6 DEVIATION-1); APD6/DRAMP arms untouched",
        "catalog_size": len(corpus), "changed": changed, "results": rows}, indent=1))
    print(f"catalog={len(corpus)} changed={len(changed)}", flush=True)
    for c in changed:
        print(json.dumps(c), flush=True)

if __name__ == "__main__":
    main()
