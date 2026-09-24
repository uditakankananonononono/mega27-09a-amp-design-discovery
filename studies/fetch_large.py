"""Scale-up fetch: full UniProt KW-0929 AMP set (10-150 aa) + reviewed non-AMP pool."""
import json, time, urllib.parse, urllib.request, hashlib, os, sys

RAW = os.path.join(os.path.dirname(__file__), "..", "data", "raw")
os.makedirs(RAW, exist_ok=True)
BASE = "https://rest.uniprot.org/uniprotkb/stream"

def stream(query):
    url = BASE + "?query=" + urllib.parse.quote(query) + "&format=fasta"
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "mega27-09a-amp/1.0"})
            with urllib.request.urlopen(req, timeout=180) as r:
                return url, r.read().decode("utf-8", "replace")
        except Exception:
            if attempt == 3: raise
            time.sleep(2 ** attempt * 2)

manifest = []
for name, q in [
    ("uniprot_kw0929_all_10-150.fasta", "keyword:KW-0929 AND length:[10 TO 150]"),
    ("uniprot_nonamp_pool_10-150.fasta", "reviewed:true NOT keyword:KW-0929 NOT keyword:KW-0800 NOT keyword:KW-0862 AND length:[10 TO 150]"),
]:
    url, text = stream(q)
    p = os.path.join(RAW, name)
    open(p, "w").write(text)
    n = text.count("\n>")
    manifest.append({"file": name, "url": url, "records": n, "bytes": len(text),
                     "sha256": hashlib.sha256(text.encode()).hexdigest()[:16],
                     "fetched_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())})
    print(name, n, "records", flush=True)
json.dump(manifest, open(os.path.join(RAW, "large_provenance.json"), "w"), indent=2)
print("DONE")
