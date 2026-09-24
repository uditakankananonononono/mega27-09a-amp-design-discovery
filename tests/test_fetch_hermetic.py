"""Hermetic checks on the fetch module (parsing only, no network)."""
from pepdesign.data.fetch import _parse_fasta


def test_parse_fasta_multiline():
    text = ">a\nACD\nEFG\n>b\nKKK\n"
    seqs = list(_parse_fasta(text))
    assert seqs == ["ACDEFG", "KKK"]


def test_parse_fasta_empty():
    assert list(_parse_fasta("")) == []
