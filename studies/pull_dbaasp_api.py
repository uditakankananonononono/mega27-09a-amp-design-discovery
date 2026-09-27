"""Q3 DBAASP full-catalog pull (addendum 6 + DEVIATION-1, locked 2026-09-27).

DEVIATION-1 (disclosed): the in-repo dbaasp_all.fasta is DEGENERATE - 546
entries but only 12 unique sequences (original novelty-ref pull bug, headers
carry repeated parent ids). Universe redefined BEFORE any label fetch: all
canonical 20-aa 5-150 aa monomer sequences from the live DBAASP catalog
(/peptides?offset&limit=100, 25,542 records). Pass 1 indexes sequences (no
activity data in list payloads). Pass 2 fetches every catalog record's full
JSON (with targetActivities) into data/raw/dbaasp_api/ - pull is blind: no
label inspection. Also dumps the corrected reference FASTA
data/raw/novelty_refs/dbaasp_catalog_2026-09-27.fasta (+ sha256) to replace
the degenerate novelty reference.
"""
import hashlib, json, pathlib, time, urllib.request

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "data/raw/dbaasp_api"
CAT = ROOT / "data/raw/novelty_refs/dbaasp_catalog_2026-09-27.fasta"
MAN = ROOT / "docs/DBAASP_API_MANIFEST_2026-09-27.md"
UA = {"User-Agent": "research-pull/1.0"}
AA20 = set("ACDEFGHIKLMNPQRSTVWY")

def get(url):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read())

INDEX = OUT / "_catalog_index.json"
OUT.mkdir(parents=True, exist_ok=True)

if INDEX.exists():
    idx = json.loads(INDEX.read_text())
    seq2id = {s: i for s, i in idx["seq2id"]}
    all_ids = idx["all_ids"]
    total = idx["total"]
    print(f"loaded index: {total} records, {len(seq2id)} sequences", flush=True)
else:
    seq2id, all_ids, offset, total = {}, [], 0, None
    while True:
        d = get(f"https://dbaasp.org/peptides?offset={offset}&limit=100")
        total = d["totalCount"]
        for rec in d["data"]:
            all_ids.append(rec["id"])
            s = rec.get("sequence")
            if s:
                seq2id.setdefault(s, rec["id"])
            for m in rec.get("monomers") or []:
                ms = m.get("sequence")
                if ms:
                    seq2id.setdefault(ms, rec["id"])
        offset += 100
        if offset % 2000 == 0:
            print(f"catalog {offset}/{total} seqs={len(seq2id)}", flush=True)
        if offset >= total:
            break
        time.sleep(0.15)
    INDEX.write_text(json.dumps({"total": total, "all_ids": all_ids,
                                 "seq2id": list(seq2id.items())}))
    print(f"catalog done: {total} records, {len(seq2id)} sequences", flush=True)

# corrected reference FASTA: canonical 5-150 aa monomer sequences
canon = sorted(s for s in seq2id if 5 <= len(s) <= 150 and set(s) <= AA20)
CAT.write_text("".join(f">DBAASP_{seq2id[s]}_{i}\n{s}\n" for i, s in enumerate(canon)))
h = hashlib.sha256(CAT.read_bytes()).hexdigest()
print(f"catalog FASTA: {len(canon)} canonical sequences sha256 {h}", flush=True)

done, failed = 0, []
for i, pid in enumerate(all_ids):
    p = OUT / f"{pid}.json"
    if p.exists() and p.stat().st_size > 100:
        done += 1
        continue
    try:
        d = get(f"https://dbaasp.org/peptides/{pid}")
        p.write_text(json.dumps(d))
        done += 1
    except Exception as e:
        failed.append((pid, str(e)))
        time.sleep(2)
    if i % 500 == 0:
        print(f"full {i}/{len(all_ids)} done={done} failed={len(failed)}", flush=True)
    time.sleep(0.4)

MAN.write_text("# DBAASP API pull manifest (addendum 6, DEVIATION-1)\n"
               "pull date 2026-09-27; universe = full DBAASP catalog "
               f"({total} records); corrected reference FASTA "
               f"dbaasp_catalog_2026-09-27.fasta {len(canon)} canonical "
               f"sequences sha256 {h}\n"
               f"full records fetched: {done}; failed: {[f[0] for f in failed]}\n"
               "pull blind: labels stored, not inspected\n")
print(f"DONE fetched={done} failed={len(failed)}", flush=True)
