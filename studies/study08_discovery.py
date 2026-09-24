"""Study 8 - AMP discovery screen over reviewed uncharacterized proteins.

Loads the study07 ensemble checkpoints (PepCNN + PepGNN trained on the
homology-aware large split) and scores reviewed UniProt 'uncharacterized
protein' entries (10-60 aa, KW-0929 excluded at fetch). Reports the top
candidates by ensemble score with novelty vs the AMP training corpus,
length, and UniProt accession. This is an in-silico prioritization screen;
candidates are hypotheses for downstream validation, not confirmed AMPs.
"""
from __future__ import annotations
import csv, json, pathlib

import numpy as np
import torch

from pepdesign.encoding import combined_features
from pepdesign.models.cnn import PepCNN
from pepdesign.models.gnn import PepGNN
from pepdesign.models.ensemble import predict_scores, predict_scores_shared_adj
from pepdesign.design import novelty

ROOT = pathlib.Path(__file__).resolve().parent.parent
TOP_N = 25


def parse_fasta(path):
    recs, header, seq = [], None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if header is not None:
                recs.append((header, "".join(seq)))
            header, seq = line[1:], []
        elif line:
            seq.append(line)
    if header is not None:
        recs.append((header, "".join(seq)))
    return recs


def cnn_inputs(seqs, max_len):
    return (torch.from_numpy(np.stack([combined_features(s, max_len) for s in seqs])),)


def accession(header):
    parts = header.split("|")
    return parts[1] if len(parts) > 1 else header.split()[0]


def main():
    cnn_ckpt = torch.load(ROOT / "results/study07_cnn.pt", weights_only=False)
    gnn_ckpt = torch.load(ROOT / "results/study07_gnn.pt", weights_only=False)
    max_len = cnn_ckpt["max_len"]
    cnn = PepCNN(); cnn.load_state_dict(cnn_ckpt["state_dict"]); cnn.eval()
    gnn = PepGNN(); gnn.load_state_dict(gnn_ckpt["state_dict"]); gnn.eval()

    recs = parse_fasta(ROOT / "data/raw/uniprot_uncharacterized_10-60.fasta")
    seqs = [s for _, s in recs]
    print(f"screening {len(seqs)} uncharacterized proteins (10-60 aa)", flush=True)

    s1 = predict_scores(cnn, cnn_inputs, seqs, max_len)
    s2 = predict_scores_shared_adj(gnn, seqs, max_len)
    ens = 0.5 * (s1 + s2)

    # novelty reference: AMP positives from the large train split
    corpus = []
    with open(ROOT / "data/processed_large/train.csv") as f:
        for row in csv.DictReader(f):
            if row["label"] == "1":
                corpus.append(row["sequence"])
            if len(corpus) >= 3000:
                break

    order = np.argsort(-ens)
    hits = []
    for i in order[:TOP_N]:
        h, s = recs[i]
        hits.append({"accession": accession(h), "header": h[:120], "sequence": s,
                     "length": len(s), "ensemble_score": float(ens[i]),
                     "cnn_score": float(s1[i]), "gnn_score": float(s2[i]),
                     "novelty_vs_amp_corpus": novelty(s, corpus)})
    for hit in hits[:10]:
        print(f"  {hit['accession']} score={hit['ensemble_score']:.3f} "
              f"novelty={hit['novelty_vs_amp_corpus']:.2f} len={hit['length']} {hit['sequence']}",
              flush=True)

    out = ROOT / "results/study08_discovery.json"
    out.write_text(json.dumps({"screened": len(seqs),
                               "max_len": max_len, "top": hits}, indent=2))
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
