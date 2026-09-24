"""Tiny bundled fixtures for hermetic tests: real published peptide sequences,
kept small so CI never touches the network."""
# Real AMPs from APD3/UniProt (short canonical sequences, positive class)
AMP_FIXTURE = [
    "GIGKFLHSAKKFGKAFVGEIMNS",   # magainin 2 (Xenopus laevis, UniProt P11006)
    "GLFDIVKKVVGALGSL",          # aurein 1.2 (UniProt P69032)
    "KWKLFKKIEKVGQNIRDGIIKAGPAVAVVGQATQIAK",  # cecropin A (P50723)
    "FFGHLFKLATKIIPSLFQ",        # temporin-derived AMP
    "GLLSGVLGVGKKIVCGLSGLC",     # dermaseptin family member
    "ILGKLLSTAWGLLSKL",          # synthetic analog of maculatin
    "GMASKAGAIAGKIAKVALKAL",     # PGLa (P11007)
    "FLPIIGKLLSGLSGLL",          # brevinin fragment
]
# Real non-AMP peptides (negative class): ordinary protein fragments
NON_AMP_FIXTURE = [
    "MKTAYIAKQRQISFVKSHFSRQ",
    "GAVLIPFYWSTCMNQDEKRHGAV",
    "SSDTPEKTLSAQVQAAVDRLG",
    "AVPATAAQGDVAQRLSQLETG",
    "KQITELAHFVRGTDKLSPQNV",
    "DVELKNGTLEAQVSAVFGHPQR",
    "TADQGRELPWSYVGNVKDAFE",
    "PLSEQAKLVGDYNATVGRWLD",
]
