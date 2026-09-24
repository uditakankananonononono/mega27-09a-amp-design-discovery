"""Hermetic tests for the pepdesign screen CLI."""
import json
import numpy as np
import torch

from pepdesign.cli import read_fasta, screen, main
from pepdesign.models.cnn import PepCNN
from pepdesign.models.gnn import PepGNN


def _write_ckpts(tmp_path):
    cnn, gnn = PepCNN(), PepGNN()
    torch.save({"state_dict": cnn.state_dict(), "max_len": 150}, tmp_path / "c.pt")
    torch.save({"state_dict": gnn.state_dict(), "max_len": 150}, tmp_path / "g.pt")
    return str(tmp_path / "c.pt"), str(tmp_path / "g.pt")


def test_read_fasta(tmp_path):
    p = tmp_path / "x.fa"
    p.write_text(">a desc\nACDEFG\nHIK\n>b\nLLLL\n")
    assert read_fasta(str(p)) == [("a", "ACDEFGHIK"), ("b", "LLLL")]


def test_screen_ranks_and_audits(tmp_path):
    c, g = _write_ckpts(tmp_path)
    ref = [("r1", "MKKLALILFMGTLVSFYADAGRKPCSGSKGGISHCTAGGKFVCNDGSISASKKTCTN")]
    cands = [("same", ref[0][1]), ("diff", "DEFWYKDEFWYKDEFWYK")]
    rows = screen(cands, c, g, ref)
    assert len(rows) == 2
    same = [r for r in rows if r["name"] == "same"][0]
    assert same["max_3mer_jaccard_vs_reference"] == 1.0
    assert same["novelty_3mer"] < 0.05
    assert all(0.0 <= r["ensemble_score"] <= 1.0 for r in rows)


def test_main_writes_json(tmp_path):
    c, g = _write_ckpts(tmp_path)
    fa = tmp_path / "cands.fa"
    fa.write_text(">x\nMKKLALILFMGTLVSFYADAGRKPCSGSKGG\n>y\nDEFWYKDEFWYKDEFWYK\n")
    ref = tmp_path / "ref.fa"
    ref.write_text(">r\nMKKLALILFMGTLVSFYADAGRKPCSGSKGGISHCTAGGKFVCNDGSISASKKTCTN\n")
    out = tmp_path / "out.json"
    rc = main([str(fa), "--cnn", c, "--gnn", g, "--reference", str(ref),
               "--out", str(out)])
    assert rc == 0
    d = json.loads(out.read_text())
    assert d["screened"] == 2 and len(d["top"]) == 2
