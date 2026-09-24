"""Study 18 - compositional audit of the large training set by label and review status.
Mechanistic explanation of the review-status confound: compares net charge,
hydropathy, length and aromatic fraction across the four label x review-status
groups. Review status is read from the FASTA header prefix (sp| reviewed, tr|
unreviewed), matching study 17."""
import csv, json, pathlib
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]
import sys; sys.path.insert(0, str(ROOT))
from pepdesign.encoding import HYDROPHOBICITY, CHARGE, AROMATIC

groups = {}
for split in ('train', 'val', 'test'):
    for r in csv.DictReader(open(ROOT/f'data/processed_large/{split}.csv')):
        label = 'pos' if r['label'] == '1' else 'neg'
        rev = 'reviewed' if r['header'].startswith('sp|') else 'unreviewed'
        groups.setdefault((label, rev), []).append(r['sequence'])

out = {}
for (label, rev), seqs in sorted(groups.items()):
    n = len(seqs)
    charge = [sum(CHARGE.get(a, 0.0) for a in s) for s in seqs]
    hydro = [np.mean([HYDROPHOBICITY[a] for a in s]) if s else 0.0 for s in seqs]
    arom = [np.mean([1.0 if a in AROMATIC else 0.0 for a in s]) if s else 0.0 for s in seqs]
    length = [len(s) for s in seqs]
    out[f'{label}_{rev}'] = {
        'n': n,
        'net_charge_mean_sd': [float(np.mean(charge)), float(np.std(charge))],
        'kd_mean_sd': [float(np.mean(hydro)), float(np.std(hydro))],
        'aromatic_frac_mean_sd': [float(np.mean(arom)), float(np.std(arom))],
        'length_mean_sd': [float(np.mean(length)), float(np.std(length))],
    }
json.dump(out, open(ROOT/'results/study18_compositional_audit.json', 'w'), indent=1)
for k, v in out.items():
    print(k, v['n'], 'charge', round(v['net_charge_mean_sd'][0],2), 'kd', round(v['kd_mean_sd'][0],2),
          'len', round(v['length_mean_sd'][0],1))
