"""Fetch reviewed UniProt 'uncharacterized protein' entries (10-60 aa) for the
AMP discovery screen. Live call; cached to data/raw/ with provenance."""
import json, time, urllib.parse, urllib.request, hashlib, os

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

q = ('reviewed:true AND length:[10 TO 60] AND '
     '(protein_name:"uncharacterized protein" OR protein_name:"hypothetical protein") '
     'NOT keyword:KW-0929')
url, text = stream(q)
p = os.path.join(RAW, "uniprot_uncharacterized_10-60.fasta")
open(p, "w").write(text)
n = text.count("\n>") + (1 if text.startswith(">") else 0)
manifest = [{"file": "uniprot_uncharacterized_10-60.fasta", "url": url,
             "records": n, "bytes": len(text),
             "sha256": hashlib.sha256(text.encode()).hexdigest()[:16],
             "fetched_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}]
json.dump(manifest, open(os.path.join(RAW, "uncharacterized_provenance.json"), "w"), indent=2)
print("records:", n, "bytes:", len(text))
