# Pre-registration addendum - 2026-09-26 (rules 1-7 revival)
Locked BEFORE any study-19 outcome is scored. Original registrations unchanged; this
addendum governs only the new work below. Amendments after outcomes will be logged as
new addenda, never edits.

## Context
Study 14 (full-record Feb2020 audit) placed the frozen ensemble at mean MCC
0.8076 / ACC 0.9033 / AUROC 0.9636 vs the published AMP Scanner Vr.2 Feb2020
10-fold CV row (SENS 0.906, SPEC 0.891, ACC 0.899, MCC 0.799, AUROC 0.962) -
numerically above but inside noise. User rule (2026-09-26 4:11:49): each project
must BEAT its benchmark - improve until it does. Study 19 is the locked
improvement attempt.

## Study 19 - locked protocol
- Identical folds to study 14: StratifiedKFold(10, shuffle=True, random_state=42),
  all 4,042 original Feb2020 records, X treated as unknown.
- Identical base learners to study 14 (2x PepCNN seeds 7/27, 1x PepGNN seed 11,
  RF-rich seed 42) PLUS one new branch: histogram gradient boosting on
  rich_features_v2 (extended physicochemical set: adds aromaticity, aliphatic
  fraction, estimated pI, max-window mean KD, grouped-AA composition,
  hydrophobic moment at 100 and 160 degrees).
- Locked a-priori ensemble weights, chosen BEFORE outcomes:
  CNN7 0.20, CNN27 0.20, GNN 0.15, RF 0.20, HGB 0.25. Threshold 0.5. No tuning.

## Locked acceptance gates (beat = ALL of G1-G3)
- G1: mean MCC across the 10 folds > 0.799 AND lower bound of the one-sided 95%
  t-interval across folds > 0.799.
- G2: mean ACC > 0.899.
- G3: mean AUROC > 0.962.
Any failed gate = documented negative, then rule-6 redirection (ChatGPT pivot
options) and a new locked addendum. No re-weighting after outcomes.

## Discovery arm (unchanged scope)
After gates, the improved ensemble re-screens the study-12 candidate set and the
uncharacterized UniProt pool; named candidates get evidence cards (score,
novelty, signal-peptide ablation, hemolysis/toxicity flags). A candidate is a
discovery only under the original evidence-card rules; zero candidates =
documented negative + pivot.

## Judge rounds (user rule: minimum 10, weaknesses + improvements)
ChatGPT adversarial rounds logged verbatim in JUDGE_ROUNDS.md with conversation
URLs; each round = weakness found, fix applied, code/data evidence of fix.
ChatGPT output is untrusted advice: factual claims are independently verified
before adoption.
