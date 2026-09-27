# Preregistration Addendum 8 (2026-09-27, locked BEFORE any Q8/Q9/Q10 scoring)

Scope: execution locks for Q8, Q9, Q10 from addendum 4. None of these arms has
been computed or inspected before this lock. Q11 (external cohort) is scoped
separately and will get its own lock once the cohort data are fixed.

## Q8 (=A14/W12): Alignment-based novelty metric for discovery candidates

- Engine: MMseqs2 easy-search (static binary already used in study23/27),
  locked params: -s 7.5, --max-seqs 10000, e-value cutoff 1e-3, no min-seq-id
  filter (alignment evidence is reported, not thresholded away).
- Query set: the 10 named InstiAMP-21 candidates + 9 prior study12 candidates
  from results/study21_rescreen.json (post-rerun version, 132df1a).
- Targets: the three pinned novelty corpora of A2.2 - APD6 (3,306), DRAMP
  (12,816), DBAASP corrected catalog (16,239, sha 735a966c).
- Reported per candidate per corpus: best-hit pident, e-value, query coverage,
  target id/sequence. Locked novelty flag: NOVEL if no hit with pident >= 30%
  AND query coverage >= 0.5 at e <= 1e-3 in any corpus.
- Comparator of record: the 3-mer Jaccard tier audit (study21 rerun). Tier/flag
  agreement is reported; any disagreement is reported as-is (neither overrides
  the other; both are pre-registered lenses).
- Output: results/study29_q8_alignment_novelty.json.

## Q9 (=A11): Negative controls through the discovery screen

- Arms: (a) scrambled - each of the 10 named InstiAMP-21 sequences shuffled
  10x with rng(11), preserving length and composition; (b) known non-AMP -
  50 reviewed non-AMP proteins drawn from the Q1 reviewed-negative pool
  (data/raw/q1_matched_negatives.fasta), rng(11) sample, lengths 10-150 aa.
- Screen: the locked study21 rerun ensemble scoring path (same weights,
  same branch recipe, retrained identically or branch-dumped predictions
  reused where available; mechanics disclosed in the result commit).
- Locked expectation, not a success threshold: scrambles should score low;
  known non-AMPs should score low. Report full p-distributions per arm
  (median, IQR, max) and every p >= 0.5 case individually, as-is. High-scoring
  controls are reported, never dropped.
- Output: results/study30_q9_negative_controls.json.

## Q10 (=A12): Encoding ablation (one-hot vs physicochemical channels)

- Dataset/protocol: feb2020 2021+2021, identical study19 rng(11) random
  10-fold splits, identical training recipe; only the input encoding varies.
- Arms (locked): (a) full combined encoding (comparator of record = study19
  numbers), (b) one-hot-only, (c) physicochemical-channels-only. Applied to
  the CNN branches (the encoding-sensitive part); RF/HGB branches unchanged
  and reported separately as context.
- Metrics: OOF MCC primary, AUROC secondary. Delta vs full encoding reported
  as-is; no success threshold; a drop is the expected citable outcome.
- Output: results/study31_q10_encoding_ablation.json.

## Rules (carried)
Negatives are reported, never terminal; deviations disclosed in-addendum before
use; each item commits with its own hash.
