# PREREG ADDENDUM 4 - 2026-09-27: LOCK of the user-verdict amendment queue (20 weaknesses + 20 additions)

Locked BEFORE any amendment execution, per the verdict protocol. Source verdict:
docs/JUDGE_VERDICT_USER_2026-09-27.md (wamid.HBgMOTE4MTM0MDk4NTcxFQIAEhgWM0VCMEQ0OUM1NzUyMTgxOTJDQzE2NwA=,
author=user, 11:26:52 IST). No outcomes for any item below have been seen at lock time
except where marked PARTIAL (pre-planned work already locked under addenda 1-3).

## Overlap map (honest)
- W7 random-split homology leakage / A3 MMseqs2 clustering: study21 hard family split
  (addendum-3, locked 04:00 AM, t0.40 DONE STRONG, t0.35 DONE HARD-pass, t0.45 running)
  is PARTIAL - it strengthens homology control with 3-mer Jaccard single-linkage but is
  NOT alignment-based MMseqs2. The MMseqs2 arm is NEW work (queue item Q2).
- W10 10-mer weak control: same PARTIAL via study21; MMseqs2 arm covers it fully.
- W5 signal-peptide ablation: already executed pre-verdict (existing result); no new work.
- W12 novelty metric: acknowledged limitation in study20/21 docs; A14 (BLAST/HMMER
  alignment-based novelty) is NEW work (queue Q8).
- Everything else: NEW.

## Amendment queue (execution order; compute side, this agent)
- Q1 (=A1, TOP PRIORITY): Reviewed-only rebuild. Fetch UniProt reviewed (Swiss-Prot)
  KW-0929 AMP positives; match negatives on review status (reviewed only) AND organism
  distribution; rebuild splits with the same study19 ensemble and rng(11) 10-fold
  protocol; retrain; re-evaluate. LOCKED REPORTING: the clean-rebuild AUROC/MCC is
  reported against the confounded 0.921/0.926 with the delta as the headline; NO success
  threshold is pre-stated - a performance DROP is the expected, citable outcome
  (the confound audit, A17, is the primary contribution). Any direction of result is
  reported as-is.
- Q2 (=A3/W10): MMseqs2 alignment-based clustering arm (static binary, locked params:
  min seq id 0.3, cov 0.8, single-linkage) replacing the 10-mer/Jaccard control;
  same folds protocol. PARTIAL-overlap noted above.
- Q3 (=A2): DBAASP potency regression (continuous MIC target from the sha256-verified
  DBAASP copy already in-repo); lock target transform (log MIC), model, and metrics
  (Spearman, MAE) before training.
- Q4 (=A7): Per-family performance report on existing predictions (anionic,
  disulfide-rich, etc. by motif/annotation) - cheap, existing data.
- Q5 (=A8): Calibration (Platt + isotonic) on held-out folds; report Brier + slope.
- Q6 (=A9): Ensemble weights optimized on validation folds vs simple average; locked
  comparison, paired folds.
- Q7 (=A13): Length-confound test: length-matched reviewed-only subset retrain.
- Q8 (=A14/W12): Alignment-based novelty metric (BLAST/HMMER vs public AMP sets) for
  discovery candidates.
- Q9 (=A11): Negative controls: scrambled sequences + known non-AMP proteins scored
  by the discovery screen.
- Q10 (=A12): Encoding ablation (one-hot vs physicochemical channels).
- Q11 (=A10): External held-out cohort not from UniProt/APD6 (DBAASP/recent set).
- Protocol/proposal sections ONLY (no claims): A5 wet-lab collaboration protocol for
  top 3-5 candidates; A15 web-tool proposal; A6 discovery-screen pre-registration
  template adopted going forward (the 852-protein screen is relabeled post-hoc
  exploratory, per W20).
- Paper-side (builder, needs my outputs): A16 12-slide condensation; A17 lead-with-audit
  reframe (AUROC 0.921/0.926 relabeled confounded measurements); A18 AMPlify/Macrel
  comparison; A19 benchmark-standard proposal (review-status cross-tabs, length
  matching, cluster splits); A20 practitioner decision tree.
- W14 compute budget (2 CPU/2GB): disclosed as a limitation; no claim changes.

## Rules
Negatives are reported, never terminal; any deviation from this queue is disclosed
in-addendum before use; each item commits with its own hash.
