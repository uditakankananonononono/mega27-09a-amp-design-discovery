"""Regression tests for the audited benchmark and discovery methods."""
import importlib.util
from pathlib import Path
import numpy as np
from pepdesign.design import novelty

ROOT = Path(__file__).resolve().parents[1]

def load(name):
    spec=importlib.util.spec_from_file_location(name, ROOT/'studies'/f'{name}.py')
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module

def test_net_charge_counts_every_residue():
    m=load('study12_candidates')
    assert m.net_charge('KKKK') > m.net_charge('K')+2.9
    assert m.net_charge('DDDD') < m.net_charge('D')-2.9

def test_novelty_is_pairwise_not_unseen_fraction():
    seq='AAAAACCCCC'
    refs=['AAAAA', 'CCCCC']
    assert novelty(seq, refs) > 0
    assert np.isclose(novelty('AAAAA', refs), 0.0)

def test_full_record_fasta_loader_retains_unknown_and_long(tmp_path):
    m=load('study14_feb2020_full')
    fa=tmp_path/'x.fa'; fa.write_text('>long\n'+'A'*183+'\n>unknown\nAXA\n')
    assert m.load_fasta_seqs(fa)==['A'*183,'AXA']
