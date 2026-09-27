# DBAASP novelty-arm re-audit (2026-09-27) - study21 rescreen correction

Cause: data/raw/novelty_refs/dbaasp_all.fasta was degenerate (546 entries, 12
unique sequences; original pull bug). study21's tiered novelty audit used it as
the DBAASP arm, understating DBAASP similarity for all screened candidates.

Fix: corrected reference data/raw/novelty_refs/dbaasp_catalog_2026-09-27.fasta
- 16,239 canonical 5-150 aa monomer sequences from the live full DBAASP catalog
(25,542 records), sha256 735a966cfe302edfd75098a68dc7ec713dc4291ea0475390c15e7330091e00f0.

Re-audit (studies/study21_dbaasp_reaudit.py, results/study21_dbaasp_reaudit.json):
- All 19 stored candidates (10 named + 9 prior rescored) recomputed against the
  corrected corpus. APD6/DRAMP arms untouched.
- Every DBAASP arm value changed (as expected), rising from ~0-0.038 to
  0.053-0.2. Worst new max-3mer-Jaccard vs DBAASP: 0.2 (InstiAMP-21-8).
- TIER FLIPS: 0. All 19 candidates remain T1 (overall max jaccard < 0.35).
  The study21 novelty tiers and headline claims STAND.
- Residual: pool-level top40_tier_counts (T1T2: 37, T3: 0) were computed over
  the 852-pool screen; the top-40 set is not stored in the JSON, so pool-level
  re-derivation requires a rescreen rerun. Given the worst candidate-level
  value is 0.2 and the tier boundary is 0.35, flips are unlikely but the
  rerun is queued for completeness and will be reported either way.
