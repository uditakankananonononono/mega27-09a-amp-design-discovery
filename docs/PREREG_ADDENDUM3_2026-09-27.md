# Preregistration Addendum 3 (09a) - locked 2026-09-27, BEFORE any harder-split run
# Origin: study20's own honest caveat (committed in results/study20_cluster_holdout.json
# and reported to program lead): at 3-mer Jaccard >= 0.6 the holdout produced 4,001
# near-singleton clusters, making the "hard" test barely harder than the random
# split (mean MCC 0.8369 vs random-split 0.8328). This addendum locks the REAL
# hard family split. Not a judge round; no outcomes seen for any threshold below.

## A3.1 Harder family split (locked parameters)
- Same single-linkage 3-mer Jaccard clustering as study20, threshold LOCKED at
  0.40 (inside the pre-stated 0.35-0.45 harder band; 0.35 and 0.45 run as
  secondary sensitivity arms with identical protocol, clearly labeled).
- 10-fold cluster-held-out CV, whole clusters per fold, no cluster crosses
  train/test; fold assignment seeded rng(11) exactly as study20.
- Identical study19 ensemble (RF .20/HGB .25 weights, threshold 0.5), retrained
  per fold. No architecture or data changes.
- Diagnostics reported: n_clusters, largest-cluster fraction, mean cluster size.

## A3.2 Locked success criteria (no post-hoc edits)
- HARD-TEST SUCCESS: mean MCC >= 0.799 (AMP Scanner Vr.2 Feb2020 published
  10-fold CV row) on the Jaccard-0.40 holdout, with per-fold values and
  one-sided 95% lower bounds reported for all metrics.
- STRONG SUCCESS: mean MCC >= 0.8328 (study19 random-split mean).
- Degradation vs study20 is EXPECTED and is the point of the harder test; it is
  reported, not tuned away. mean MCC < 0.799 = documented negative -> rule-6
  redirection consult (judge), never silent threshold shopping.
