# Preregistration Addendum 2 (09a) - locked 2026-09-26, BEFORE any cluster-holdout
# or tiered-novelty run. Origin: ChatGPT judge round 1 (JUDGE_ROUNDS.md R1,
# conversation https://chatgpt.com/c/6ab7e420-12e0-83e8-b100-367851b30c05).
# Rule-8 fold-back: the judge's critique produced these two concrete novelty
# additions. ChatGPT output is untrusted advice; every design choice below was
# independently checked for implementability in this repo before locking.

## A2.1 Cluster-held-out generalization evaluation (NEW experiment)
Motivation (judge's sharpest hit): random-split CV can reward recognizing known
AMP families; a discovery tool must generalize to sequence families absent from
training.
Locked design:
- Clustering: greedy single-linkage over the Feb2020 dataset (2021 AMP + 2021
  decoy), edge when 3-mer Jaccard >= 0.6 (for 20-60 aa peptides this
  approximates ~70-80% pairwise identity; approximation documented in paper).
- Split: 10-fold cluster-held-out CV - folds contain whole clusters only;
  no cluster crosses train/test. Fold assignment seeded rng(11).
- Model: the identical locked study19 ensemble (CNN7 .20/CNN27 .20/GNN .15/
  RF .20/HGB .25, threshold 0.5), retrained per fold exactly as in study19.
- Report: per-fold and mean MCC/ACC/AUC with one-sided 95% t-intervals.
- Locked interpretation criterion (NOT a re-beat gate): generalization is
  claimed only if cluster-holdout mean MCC >= study19 random-split mean MCC
  (0.8328) - 0.10. A larger drop = documented negative + rule-6 redirection.

## A2.2 Tiered novelty audit for the discovery arm (TIGHTENED criterion)
Motivation: "novel" candidates that are near-duplicates of known AMPs are
rediscovery, not discovery.
Locked design:
- Every re-screen candidate is audited against APD6 2024 (in-repo) plus the
  DRAMP and DBAASP public dumps (pinned download, sha256 recorded).
- Identity proxy: 3-mer Jaccard and longest common substring vs every
  reference sequence (exact metrics already used in study12).
- Tiers: T1 max-Jaccard < 0.35 (genuinely distant), T2 0.35-0.6, T3 > 0.6
  (near-duplicate). Discovery candidates are named ONLY from T1/T2; T3 hits
  are reported as rediscovery controls.
- Evidence cards (score, novelty tiers, physicochemical profile,
  signal-peptide ablation, hemolysis flags) unchanged from addendum 1.
