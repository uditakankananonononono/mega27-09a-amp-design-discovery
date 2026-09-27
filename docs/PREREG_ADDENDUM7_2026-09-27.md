# Preregistration Addendum 7 (2026-09-27, locked BEFORE any Q7 scoring)

Scope: execution lock for Q7 (=A13) from addendum 4. No Q7 labels, metrics, or
outcomes have been computed or inspected before this lock.

## Q7 (=A13): Length-confound test - exact-length-matched reviewed-only retrain

- Data: the study22 locked Q1 dataset, unchanged:
  data/raw/q1_reviewed_positives.fasta (2,973 reviewed KW-0929 positives) +
  data/raw/q1_matched_negatives.fasta (2,973 A4.1-matched reviewed negatives).
- Matching rule (locked): exact-length 1:1 matching WITHOUT replacement. For
  each positive (in rng(11)-shuffled order), take a same-length unused negative
  (ties among negatives broken by rng(11)). Positives with no remaining
  exact-length negative are DROPPED. Report: matched pairs, dropped positives,
  unused negatives, and the matched length distribution.
- Pre-locked fallback (disclosed before use, in the A4.1 pattern): if exact
  matching yields fewer than 1,000 pairs, widen to length +/- 1 aa and report
  the fallback as used. No other widening.
- Protocol: identical study19 locked ensemble (CNN7 .20 / CNN27 .20 / GNN .15 /
  RF .20 / HGB .25, threshold 0.5), same run_fold machinery, rng(11) random
  10-fold splits (same make_splits construction as study22) on the matched
  subset only. OOF metrics computed exactly as study22.
- Metrics: MCC primary, AUROC secondary (OOF).
- Comparator of record: study22 FINAL headline (AUROC 0.9633, MCC 0.803).
  The delta vs study22 IS the length-confound estimate. No success threshold;
  a performance DROP is the expected, citable outcome; any direction is
  reported as-is.
- Falsifier: a label-shuffle arm (rng(11)) must collapse to chance
  (|MCC| < 0.1 expected); if not, the result is flagged invalid, not tuned.
- Output: results/study28_q7_length_confound.json; committed with its own hash.

## Rules (carried from addendum 4)
Negatives are reported, never terminal; deviations disclosed in-addendum before
use; each item commits with its own hash.
