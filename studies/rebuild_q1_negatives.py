"""Rebuild Q1 matched negatives - VERBATIM recovery of the original builder
(recovered from the archived build transcript after the 2026-09-27 sandbox
wipe; reproduces the pinned sha256 f446c276658653af from
results/q1_provenance.json byte-exactly; fill counts species 1503 / genus 365 /
global 1105 match the committed provenance).

Rule (results/q1_provenance.json): 1:1 match, rng(11): species -> genus ->
global length-matched (addendum-4 amendment A4.1). Positives file must be the
byte-exact q1_reviewed_positives.fasta (sha256 prefix 58f444dadbb290f4).

Run: PYTHONPATH=. python3 studies/rebuild_q1_negatives.py
"""
import re, hashlib, collections, random, bisect, pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent

def read_fasta(p):
    recs, hdr, seq = [], None, []
    for line in open(p):
        line = line.rstrip('\n')
        if line.startswith('>'):
            if hdr: recs.append((hdr, ''.join(seq)))
            hdr, seq = line, []
        else: seq.append(line)
    if hdr: recs.append((hdr, ''.join(seq)))
    return recs

def org(h):
    m = re.search(r'OS=([^=]+?)\sOX=', h)
    return m.group(1).strip() if m else 'unknown'

pos = read_fasta(ROOT/'data/raw/q1_reviewed_positives.fasta')
neg = read_fasta(ROOT/'data/raw/uniprot_nonamp_pool_10-150.fasta')
rng = random.Random(11)
neg_by_species = collections.defaultdict(list); neg_by_genus = collections.defaultdict(list)
for h, s in neg:
    o = org(h); neg_by_species[o].append((h, s)); neg_by_genus[o.split()[0]].append((h, s))
for v in neg_by_species.values(): rng.shuffle(v)
for v in neg_by_genus.values(): rng.shuffle(v)
pos_org = collections.Counter(org(h) for h, _ in pos)
matched, used, stats = [], set(), {'species': 0, 'genus': 0, 'global': 0}
for o, n in pos_org.items():
    pool = [r for r in neg_by_species.get(o, []) if r[0] not in used]
    take = pool[:n]
    for r in take: used.add(r[0])
    matched.extend(take); stats['species'] += len(take); rem = n - len(take)
    if rem > 0:
        gpool = [r for r in neg_by_genus.get(o.split()[0], []) if r[0] not in used]
        take2 = gpool[:rem]
        for r in take2: used.add(r[0])
        matched.extend(take2); stats['genus'] += len(take2); rem -= len(take2)
    if rem > 0: stats['global'] += rem
pos_lens = sorted(len(s) for _, s in pos)
rest = [r for r in neg if r[0] not in used]
rng.shuffle(rest)
rest.sort(key=lambda r: min(abs(len(r[1]) - pos_lens[min(bisect.bisect_left(pos_lens, len(r[1])), len(pos_lens)-1)]),
                            abs(len(r[1]) - pos_lens[max(bisect.bisect_left(pos_lens, len(r[1]))-1, 0)])))
matched.extend(rest[:stats['global']])
assert len(matched) == len(pos), (len(matched), len(pos))
out = ROOT/'data/raw/q1_matched_negatives.fasta'
with open(out, 'w') as f:
    for h, s in matched: f.write(h + '\n' + s + '\n')
h = hashlib.sha256(out.read_bytes()).hexdigest()
print('negatives', len(matched), 'fill', stats)
print('sha256', h[:16], 'expect f446c276658653af', 'MATCH' if h.startswith('f446c276658653af') else 'NO')
