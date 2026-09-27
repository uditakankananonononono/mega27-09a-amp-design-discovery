# Novelty-audit reference corpora (addendum 2, A2.2) - pinned
Raw files live in data/raw/novelty_refs/ (gitignored; re-fetchable).
- DRAMP natural AMPs: dramp_natural_amps.txt, 2,177,343 bytes, retrieved
  2026-09-26 from http://dramp.cpu-bioinfor.org/downloads/download.php?filename=download_data%2FDRAMP3.0_new%2Fnatural_amps.txt
  sha256 fbaebb527695785ec3e5c0d13cf6a4eb82aadfeeb5bca64e534fcfe22efd12aa
- DBAASP: API pull https://dbaasp.org/peptides (25,542 records, pass-2 pull
  COMPLETE 2026-09-27 17:38 IST, 25,543 files incl. index);
  canonical 5-150 aa monomer sequences -> dbaasp_all.fasta, replaced
  2026-09-27 with the corrected full catalog (dbaasp_catalog_2026-09-27.fasta
  content; the pre-DEVIATION-1 degenerate file is superseded, in git history),
  16,239 sequences, 1,090,067 bytes (catalog file), sha256
  735a966cfe302edfd75098a68dc7ec713dc4291ea0475390c15e7330091e00f0
- APD6 2024: data/raw/apd6_natural_2024.fasta (already pinned in-repo history).
