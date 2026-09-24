"""Live dataset fetchers. NEVER called in CI - studies call these only when run
with --live, then cache raw FASTA/CSV under data/raw/."""
from __future__ import annotations
import io
import os
import time

import requests

UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"


def fetch_uniprot_keyword_fasta(keyword_id: str, max_records: int, out_path: str,
                                min_len: int = 8, max_len: int = 60,
                                reviewed: bool = True) -> int:
    """Stream reviewed UniProt entries carrying a keyword (e.g. KW-0929
    'Antimicrobial') as FASTA, keeping length-bounded sequences."""
    query = f"keyword:{keyword_id}"
    if reviewed:
        query += " AND reviewed:true"
    params = {"query": query, "format": "fasta", "size": min(max_records, 500)}
    kept = 0
    cursor_url = None
    with open(out_path, "w") as fh:
        while kept < max_records:
            if cursor_url is None:
                resp = requests.get(UNIPROT_STREAM, params=params, timeout=60)
            else:
                resp = requests.get(cursor_url, timeout=60)
            resp.raise_for_status()
            seqs = _parse_fasta(resp.text)
            for seq in seqs:
                if min_len <= len(seq) <= max_len and set(seq) <= set("ACDEFGHIKLMNPQRSTVWY"):
                    fh.write(f">u{kept}\n{seq}\n")
                    kept += 1
                    if kept >= max_records:
                        break
            link = resp.headers.get("Link", "")
            if 'rel="next"' not in link or kept >= max_records:
                break
            cursor_url = link.split("<")[1].split(">")[0]
            time.sleep(0.3)
    return kept


def fetch_uniprot_negatives(exclude_keyword: str, max_records: int, out_path: str,
                            min_len: int = 8, max_len: int = 60) -> int:
    """Reviewed cytosolic protein fragments lacking the keyword, length-matched."""
    query = f"reviewed:true AND NOT keyword:{exclude_keyword} AND length:[{min_len} TO {max_len}]"
    params = {"query": query, "format": "fasta", "size": min(max_records, 500)}
    resp = requests.get(UNIPROT_STREAM, params=params, timeout=60)
    resp.raise_for_status()
    kept = 0
    with open(out_path, "w") as fh:
        for seq in _parse_fasta(resp.text):
            if kept >= max_records:
                break
            if min_len <= len(seq) <= max_len and set(seq) <= set("ACDEFGHIKLMNPQRSTVWY"):
                fh.write(f">n{kept}\n{seq}\n")
                kept += 1
    return kept


def _parse_fasta(text: str):
    buf = None
    for line in text.splitlines():
        if line.startswith(">"):
            if buf is not None:
                yield "".join(buf)
            buf = []
        elif buf is not None:
            buf.append(line.strip())
    if buf:
        yield "".join(buf)
