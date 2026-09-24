"""Fetch canonical AMP Scanner datasets; validate checksums and FASTA counts."""
from __future__ import annotations
import hashlib, io, json, pathlib, urllib.request, zipfile
from Bio import SeqIO
ROOT=pathlib.Path(__file__).resolve().parent.parent
VELTRI='https://raw.githubusercontent.com/dan-veltri/amp-scanner-v2/main/original-dataset/'
FILES={'AMP.tr.fa':712,'DECOY.tr.fa':712,'AMP.eval.fa':354,'DECOY.eval.fa':354,'AMP.te.fa':712,'DECOY.te.fa':712}
FEB='https://www.dveltri.com/ascan/v2/data/AMP_Scan2_Feb2020_Dataset.zip'

def fetch(url):
    with urllib.request.urlopen(url, timeout=30) as r: return r.read()

def main():
    v=ROOT/'data/raw/veltri'; v.mkdir(parents=True, exist_ok=True)
    manifest={'veltri_source':VELTRI,'feb2020_source':FEB,'files':{}}
    for name,n in FILES.items():
        b=fetch(VELTRI+name); (v/name).write_bytes(b)
        count=sum(1 for _ in SeqIO.parse(io.StringIO(b.decode()),'fasta'))
        if count != n: raise ValueError((name,count,n))
        manifest['files'][name]={'sha256':hashlib.sha256(b).hexdigest(),'records':count}
    b=fetch(FEB); z=zipfile.ZipFile(io.BytesIO(b)); f=ROOT/'data/raw/feb2020'; f.mkdir(exist_ok=True)
    for name in ('AMPS_02182020.fasta','DECOYS_02182020.fasta'):
        matches=[x for x in z.namelist() if x.endswith(name)]
        if len(matches)!=1: raise ValueError((name,matches))
        content=z.read(matches[0]); (f/name).write_bytes(content)
        n=sum(1 for _ in SeqIO.parse(io.StringIO(content.decode()),'fasta'))
        if n!=2021: raise ValueError((name,n))
        manifest['files'][name]={'sha256':hashlib.sha256(content).hexdigest(),'records':n}
    (ROOT/'results/benchmark_sources.json').write_text(json.dumps(manifest,indent=2))
    print(json.dumps(manifest,indent=2))
if __name__=='__main__': main()
