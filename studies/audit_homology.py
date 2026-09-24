"""Exact k-mer overlap audit between benchmark folds; not an alignment test."""
from __future__ import annotations
import json,pathlib
import numpy as np
from sklearn.model_selection import StratifiedKFold
from studies.study14_feb2020_full import load_fasta_seqs
ROOT=pathlib.Path(__file__).resolve().parents[1]

def kmer_set(s,k): return {s[i:i+k] for i in range(len(s)-k+1)}
def audit(seqs,labels,k=10):
    skf=StratifiedKFold(n_splits=10,shuffle=True,random_state=42)
    rows=[]
    for j,(tri,tei) in enumerate(skf.split(seqs,labels),1):
        train=set().union(*(kmer_set(seqs[i],k) for i in tri))
        test=[seqs[i] for i in tei]
        overlap=[bool(kmer_set(s,k)&train) for s in test]
        rows.append({'fold':j,'test_n':len(tei),'shared_exact_10mer_n':sum(overlap),
                     'fraction':sum(overlap)/len(overlap)})
    return rows

def main():
    raw=ROOT/'data/raw/feb2020'
    pos=load_fasta_seqs(raw/'AMPS_02182020.fasta')
    neg=load_fasta_seqs(raw/'DECOYS_02182020.fasta')
    seqs=pos+neg; y=np.array([1]*len(pos)+[0]*len(neg))
    rows=audit(seqs,y)
    out={'dataset_n':len(seqs),'k':10,'folds':rows,
         'overall_fraction':sum(r['shared_exact_10mer_n'] for r in rows)/len(seqs)}
    (ROOT/'results/study16_homology_audit.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
