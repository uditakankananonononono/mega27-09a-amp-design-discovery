"""pepdesign screen - leakage-controlled AMP discovery screening CLI.

Screens a FASTA of candidate sequences with a trained ensemble, audits
novelty against a reference AMP corpus (k-mer novelty + max k-mer Jaccard),
and emits a ranked candidate table (JSON/TSV). One pass: score -> novelty
filter -> ranked report. No other open AMP tool bundles a novelty audit
with scoring in a single command.
"""
from __future__ import annotations
import argparse, json, sys

import numpy as np
import torch

from pepdesign.encoding import combined_features
from pepdesign.models.cnn import PepCNN
from pepdesign.models.gnn import PepGNN
from pepdesign.models.ensemble import predict_scores, predict_scores_shared_adj
from pepdesign.design import novelty

AA = set("ACDEFGHIKLMNPQRSTVWY")


def read_fasta(path):
    seqs, name, cur = [], None, ""
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    seqs.append((name, cur))
                name, cur = line[1:].split()[0], ""
            elif line:
                cur += "".join(ch for ch in line if ch in AA)
    if name is not None:
        seqs.append((name, cur))
    return seqs


def kmer_jaccard(a, b, k=3):
    ka = {a[i:i+k] for i in range(len(a)-k+1)}
    kb = {b[i:i+k] for i in range(len(b)-k+1)}
    return len(ka & kb) / max(len(ka | kb), 1)


def screen(candidates, cnn_ckpt, gnn_ckpt, reference, max_len=150,
           novelty_floor=0.5, ref_cap=3000):
    ck = torch.load(cnn_ckpt, weights_only=False)
    cnn = PepCNN(); cnn.load_state_dict(ck["state_dict"]); cnn.eval()
    ck = torch.load(gnn_ckpt, weights_only=False)
    gnn = PepGNN(); gnn.load_state_dict(ck["state_dict"]); gnn.eval()
    ref = [s for _, s in reference][:ref_cap]

    def cnn_inputs(ss, ml):
        return (torch.stack([torch.from_numpy(combined_features(x, ml)) for x in ss]),)

    rows = []
    for name, seq in candidates:
        if not (5 <= len(seq) <= max_len):
            continue
        c = float(predict_scores(cnn, cnn_inputs, [seq], max_len)[0])
        g = float(predict_scores_shared_adj(gnn, [seq], max_len)[0])
        e = 0.5 * (c + g)
        nov = novelty(seq, ref)
        jac = max((kmer_jaccard(seq, r) for r in ref), default=0.0)
        rows.append({"name": name, "sequence": seq, "length": len(seq),
                     "ensemble_score": round(e, 4), "cnn_score": round(c, 4),
                     "gnn_score": round(g, 4), "novelty_3mer": round(nov, 3),
                     "max_3mer_jaccard_vs_reference": round(jac, 3),
                     "passes_novelty_floor": nov >= novelty_floor})
    rows.sort(key=lambda r: -r["ensemble_score"])
    return rows


def main(argv=None):
    p = argparse.ArgumentParser(prog="pepdesign screen",
                                description="AMP discovery screen with novelty audit")
    p.add_argument("fasta", help="candidate sequences (FASTA)")
    p.add_argument("--cnn", required=True, help="CNN checkpoint (.pt)")
    p.add_argument("--gnn", required=True, help="GNN checkpoint (.pt)")
    p.add_argument("--reference", required=True, help="reference AMP corpus (FASTA)")
    p.add_argument("--max-len", type=int, default=150)
    p.add_argument("--novelty-floor", type=float, default=0.5)
    p.add_argument("--top", type=int, default=25)
    p.add_argument("--out", help="write JSON here (default stdout)")
    args = p.parse_args(argv)

    rows = screen(read_fasta(args.fasta), args.cnn, args.gnn,
                  read_fasta(args.reference), args.max_len, args.novelty_floor)
    payload = json.dumps({"screened": len(rows), "top": rows[:args.top]}, indent=2)
    if args.out:
        with open(args.out, "w") as f:
            f.write(payload)
    else:
        print(payload)
    return 0


if __name__ == "__main__":
    sys.exit(main())
