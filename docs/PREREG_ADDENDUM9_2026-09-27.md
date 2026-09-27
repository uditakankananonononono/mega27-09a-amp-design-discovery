# Preregistration Addendum 9 (2026-09-27, locked BEFORE any Q11 scoring)

Scope: execution lock for Q11 (=A10) from addendum 4: external held-out cohort
not from UniProt/APD6. No Q11 scores have been computed or inspected.

## Q11 (=A10): External cohort - DBAASP activity-defined eval

- Rationale (locked design choice, disclosed before use): a non-UniProt,
  non-APD6 external cohort needs both external positives AND external
  negatives; public non-AMP negative sets are UniProt-derived by construction.
  DBAASP target-activity records give an activity-defined cohort that is
  external to every training source and tied to real experimental outcomes.
- Cohort (locked): from the completed blind DBAASP pull (25,542 records,
  manifest e4c1e5a), canonical 5-150 aa 20-letter monomers with at least one
  E. coli MIC measurement (the study27 primary-cohort record set):
  - external POSITIVES: median MIC <= 32 ug/ml (clearly active);
  - external NEGATIVES: all reported E. coli MICs > 256 ug/ml (reported
    inactive at assay range; uM converted per the study27 locked conversion).
  Cohort sizes reported; multimers excluded (study27 rule).
- Contamination control (locked): exact-sequence dedup vs ALL training data
  (Q1 positives, Q1 negatives, feb2020 AMPs + decoys) PLUS MMseqs2
  cluster-level removal - any sequence falling in a study23-param cluster
  (min-seq-id 0.3, cov 0.8, single-linkage) that touches a training sequence
  is dropped. Dropped counts reported.
- Model (locked): the study22 recipe ensemble trained on the FULL Q1 dataset
  (2,973+2,973; same 5 branches, seeds, epochs; a final-model train, the
  study21 pattern) - no threshold or weight changes.
- Metrics (locked): AUROC primary (ranking quality on external outcomes),
  MCC at the locked 0.5 threshold secondary; per-class score distributions
  reported. Spearman vs log10 median MIC on positives reported as context
  alongside study27. Reported as-is, any direction; a weak result bounds the
  Q1-clean headline's external validity and is itself the citable outcome.
- Output: results/study32_q11_external_cohort.json; own commit hash.

## Rules (carried)
Negatives are reported, never terminal; deviations disclosed in-addendum before
use; each item commits with its own hash.
