"""Q3 DBAASP full-record pull - CONCURRENT BATCH DRIVER (addendum 6).

Runs one bounded batch per invocation (~batch seconds), resume-capable via
files on disk. Sandboxes suspend between agent turns, so this is driven in
the foreground during active turns instead of a detached background wait.
Blind pull: records are written to disk without label inspection.
"""
import json, pathlib, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "data/raw/dbaasp_api"
UA = {"User-Agent": "research-pull/1.0"}
WORKERS = 6

idx = json.loads((OUT / "_catalog_index.json").read_text())
all_ids = idx["all_ids"]
todo = [i for i in all_ids
        if not ((OUT / f"{i}.json").exists() and (OUT / f"{i}.json").stat().st_size > 100)]
print(f"remaining: {len(todo)}/{len(all_ids)}", flush=True)
if not todo:
    print("PULL COMPLETE", flush=True); sys.exit(0)

budget = float(sys.argv[1]) if len(sys.argv) > 1 else 100.0
deadline = time.time() + budget

def fetch(pid):
    req = urllib.request.Request(f"https://dbaasp.org/peptides/{pid}", headers=UA)
    for t in range(3):
        try:
            with urllib.request.urlopen(req, timeout=40) as r:
                (OUT / f"{pid}.json").write_bytes(r.read())
            return True
        except Exception:
            time.sleep(2 * (t + 1))
    return False

done = fail = 0
with ThreadPoolExecutor(max_workers=WORKERS) as ex:
    futs = {}
    it = iter(todo)
    # prime
    for _ in range(WORKERS * 2):
        try: futs[ex.submit(fetch, next(it))] = None
        except StopIteration: break
    while futs:
        for fut in list(futs):
            if fut.done():
                done += 1 if fut.result() else 0; fail += 0 if fut.result() else 1
                del futs[fut]
                if time.time() < deadline:
                    try: futs[ex.submit(fetch, next(it))] = None
                    except StopIteration: pass
        time.sleep(0.05)
        if time.time() >= deadline:
            break
print(f"batch done: fetched~{done} failed~{fail}, remaining after batch ~{len(todo)-done-fail}", flush=True)
