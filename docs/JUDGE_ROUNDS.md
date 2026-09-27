# Judge rounds (rule 8: a round counts only when its critique produces a
# concrete novelty change folded back into the work). ChatGPT output is
# untrusted advice; factual claims are verified independently before adoption.

---

## RULE CHANGE 2026-09-27 (counted-round requirement 10 -> 1)
User, WhatsApp 10:00:07 IST (wamid.HBgMOTE4MTM0MDk4NTcxFQIAEhgWM0VCMDJCMTZGRTVEMkQwMTFBQzc4MQA=,
verified author=user): "NOT 10 ROUNDS OOF CHATGPT CHECK JUST ONE WHICH I PROVIDE OK?"
Effect: the counted ChatGPT judge requirement is now ONE round per project,
provided by the user through the courier route. History below is preserved, not
deleted. Status of this project's judge gate: 1 of 1 - REQUIREMENT MET (2026-09-27).
Her provided verdict arrived 11:26:52 IST (wamid.HBgMOTE4MTM0MDk4NTcxFQIAEhgWM0VCMEQ0OUM1NzUyMTgxOTJDQzE2NwA=,
verified author=user): 20 weaknesses + 20 additions, archived verbatim in
docs/JUDGE_VERDICT_USER_2026-09-27.md; amendment queue locked BEFORE execution
in docs/PREREG_ADDENDUM4_2026-09-27.md.

CLARIFICATION 2026-09-27 10:01:47 IST (user, WhatsApp
wamid.HBgMOTE4MTM0MDk4NTcxFQIAEhgWM0VCMDMwREI5RDQ0QUNCRDc2MTNDMwA=, verified
author=user): "EACH PROJECTS NEED ONE FROM ME TO PASS" - only a verdict she
personally provides through the courier paste route counts. Agent-initiated
ChatGPT rounds, even in her account, do NOT satisfy the gate. Round 1 below
(2026-09-26) was agent-initiated: it remains in this ledger as preserved
history/supplementary evidence and its landed novelty change stands as
science, but it does NOT count toward the judge gate. A courier prompt for
this project will be delivered to her via the parent; the gate passes when her
verdict returns (wamid provenance) and, per rule 8, its critique produces a
landed novelty change. Gemini/DeepSeek/LLM consults are supplementary only -
logged, never counted.

---

## Round 1 - 2026-09-26, novelty critique of the study19 state
Conversation: https://chatgpt.com/c/6ab7e420-12e0-83e8-b100-367851b30c05
Route: text paste into chatgpt.com (user rule 8:16 PM), one fill, clicked Send.

### Verbatim prompt
You are a novelty judge for a peer-review-grade antimicrobial peptide (AMP) prediction project. Verified current state: a locked-preregistration 10-fold cross-validation on the published AMP Scanner Vr.2 Feb-2020 dataset (2021 AMPs + 2021 decoys, identical fold split as the prior study) with an a-priori-locked ensemble (CNN7 .20/CNN27 .20/GNN .15/RF .20/HGB .25, threshold 0.5, HGB branch on an extended physicochemical feature set) achieved MCC 0.8328 (SD 0.0192, worst fold 0.8025) vs the published 0.799, ACC 0.9161 vs 0.899, AUC 0.9688 vs 0.962 - all three locked gates passed, every fold individually beating the published MCC. Next planned step: a discovery re-screen for novel high-confidence AMP candidates. Questions: (1) Is a locked-weights reproduction-plus-improvement claim genuinely novel science, or merely a reproduction - what would a skeptical reviewer say? (2) Give your single highest-value CONCRETE novelty addition for the next iteration: name the specific method, dataset, or analysis, implementable in days not months, that most raises the contribution. (3) Name one failure mode you would attack as a reviewer of the discovery re-screen plan. Be specific; no generic advice.

### Verbatim response (page capture; UI chrome retained)
Current URL: https://chatgpt.com/c/6ab7e420-12e0-83e8-b100-367851b30c05

Skip to content
Chat history
ChatGPT
New chat
Images
Library
Scheduled
Plugins
Projects
Codex
More
Recents
Novelty Review and Next Experiment
ISEF Review Analysis
Redirection Proposal
Awaiting Part Two
Critique RL Research Framing
Review computational biology project
Acknowledge Paper Parts
Acknowledge Research Draft
Acknowledge paper part one
CRISPR Evidence Chain Review
Critique Experimental Design
Scientific Proposal Critique
Evaluate Drug Synergy Manuscript
Propose Annotation Stability Method
Simulated ISEF Critique
Propose novelty analysis
Critical review and novelty analysis
Assess MoCo Causality
Assess CRISPR Claim Strength
Evaluate Organoid Validation
CRISPR Signal Validation Review
Scientific Review and Next Experiments
Adversarial ISEF Review
Evaluate Rescorer Validation
Adversarial Study Review
Bioinformatics Critique
Scientific Judge Review
Redirecting CRISPR Findings
Udita Kankana Phookan
Free
Claim offer
Free offer
Share
You are a novelty judge for a peer-review-grade antimicrobial peptide (AMP) prediction project. Verified current state: a locked-preregistration 10-fold cross-validation on the published AMP Scanner Vr.2 Feb-2020 dataset (2021 AMPs + 2021 decoys, identical fold split as the prior study) with an a-priori-locked ensemble (CNN7 .20/CNN27 .20/GNN .15/RF .20/HGB .25, threshold 0.5, HGB branch on an extended physicochemical feature set) achieved MCC 0.8328 (SD 0.0192, worst fold 0.8025) vs the published 0.799, ACC 0.9161 vs 0.899, AUC 0.9688 vs 0.962 - all three locked gates passed, every fold individually beating the published MCC. Next planned step: a discovery re-screen for novel high-confidence AMP candidates. Questions: (1) Is a locked-weights reproduction-plus-improvement claim genuinely novel science, or merely a reproduction - what would a skeptical reviewer say? (2) Give your single highest-value CONCRETE novelty addition for the next iteration: name the specific method, dataset, or analysis, implementable in days not months, that most raises the contribution. (3) Name one failure mode you would attack as a reviewer of the discovery re-screen plan. Be specific; no generic advice.
Show more
1) Is this genuinely novel science or merely reproduction-plus-improvement?

A skeptical peer reviewer would likely say:

The current result is a strong benchmark improvement, but by itself it is not yet a major scientific contribution.

The reason is that the scientific claim is currently:

“Using a locked ensemble and an expanded feature branch, we improve performance on the same AMP Scanner Vr.2 benchmark under the same CV protocol.”

That is methodological validation, not necessarily new biological knowledge.

A reviewer would separate the achievement into:

What is genuinely valuable

Pre-registration + locked evaluation is a major strength. Many ML-bio papers suffer from post-hoc tuning.

Beating the published baseline on all folds reduces the chance of a lucky split.

MCC 0.83 on a balanced AMP/decoy benchmark is a meaningful improvement.

The HGB physicochemical feature branch suggests that interpretable sequence chemistry contributes beyond deep architectures.

What a skeptical reviewer would criticize

They might write:

“The authors demonstrate improved discrimination on a previously curated benchmark, but the work does not establish generalization beyond the dataset construction process or reveal new antimicrobial biology.”

The main issue is benchmark saturation. AMP Scanner Vr.2 is not a blind biological universe; it is a historical labeled dataset. Improvements can arise from:

better feature engineering,

hidden dataset artifacts,

sequence-family leakage,

decoy construction biases,

organism/source differences.

So the current project is closer to:

Strong computational methods paper / high-quality validation study

rather than:

Discovery science paper

The discovery re-screen is where the novelty opportunity begins.

2) Highest-value concrete novelty addition

My highest-value addition would be:

Add a cluster-aware external generalization test using AMP family holdout + prospective candidate ranking

Specifically:

Method

Perform sequence identity clustering before evaluation using MMseqs2/CD-HIT, then create a cluster-held-out test set.

Example:

Cluster all AMP Scanner sequences at:

70% identity

50% identity

Train/test split by clusters, not random sequences.

Evaluate your locked model without retraining or with a pre-specified retraining protocol.

Compare:

Evaluation	Expected difficulty
Random 10-fold CV	easier
70% cluster holdout	harder
50% cluster holdout	hardest
Why this is the highest-value addition

Because the biggest reviewer attack is:

“Your model may recognize known AMP families rather than discover new antimicrobial peptides.”

A cluster holdout directly attacks that.

A model that maintains strong MCC/AUC on family-separated data is much closer to a discovery tool.

Implementability

This is a days-level addition:

MMseqs2/CD-HIT clustering: hours

split generation: hours

locked inference: hours

analysis/plots: 1–2 days

No wet lab required.

Bonus (very high impact)

Then use your model to screen a large unlabeled peptide universe:

Possible sources:

UniProt reviewed/unreviewed peptides

metagenomic peptide databases

antimicrobial peptide repositories without labels

Then report:

number of candidates

diversity after clustering

novelty relative to training AMP families

confidence distribution

physicochemical profile

The key is:

Do not just find high-scoring peptides. Find high-scoring peptides that are far from known AMPs.

A reviewer cares about:

“Did you discover unseen AMP-like sequence space?”

not:

“Did your classifier rank known AMP-like sequences highly?”

3) One failure mode I would attack as a reviewer
Failure mode: discovery candidates are actually database-near duplicates of known AMPs

This is the most dangerous weakness.

Imagine the pipeline:

Train model.

Scan UniProt.

Select top 100 peptides.

Call them “novel AMP candidates.”

A reviewer asks:

“Novel according to what definition?”

A peptide could be:

absent from AMP Scanner,

but 95% identical to a known AMP in another database,

or from the same AMP family,

or generated from the same organism lineage.

Then “discovery” collapses into rediscovery.

Specific attack I would perform

For every candidate:

Run:

MMseqs2 against:

DBAASP

DRAMP

APD3

UniProt AMP annotations

Report:

closest known AMP identity

alignment coverage

E-value

family membership

Then define novelty tiers:

Tier 1

<50% identity to known AMP families

Tier 2

50–80%

Tier 3

80%

Your strongest candidates should come from Tier 1/2.

Reviewer verdict today

If I reviewed the current manuscript:

Benchmark result: strong
Methodological rigor: unusually good
Biological novelty: currently limited
Discovery potential: high, but unproven

The single change most likely to move it from “better classifier” to “research contribution” is:

Cluster-held-out generalization + novelty-filtered discovery screening.

That converts the central claim from:

“We classify AMP-like sequences better”

into:

“We identify antimicrobial peptide candidates in sequence space that is genuinely separated from known AMP families.”

Do you like this personality?
ChatGPT is AI and can make mistakes. Check important info.

Think

Iframe [frame_id=EF7DBFED6E0C21E017F545DBB63F47FB parent_frame_id=071841B7B9B3C874094C898F4FC722BC origin=https://chatgpt.com depth=1]

### Critique (as used)
1. A locked reproduction-plus-improvement is strong rigor but limited
   biological novelty on its own; reviewers will attack family recognition
   vs genuine discovery.
2. Highest-value addition: cluster-held-out generalization eval +
   novelty-tiered discovery screening.
3. Failure mode: candidates that are near-duplicates of known AMPs.

### Novelty change folded back (the reason this round counts)
- docs/PREREG_ADDENDUM2_2026-09-26.md locked BEFORE any run:
  A2.1 cluster-held-out 10-fold CV of the locked study19 ensemble with a
  pre-set interpretation criterion; A2.2 tiered novelty audit (T1/T2/T3
  vs APD6 + DRAMP + DBAASP) binding on every discovery candidate.
- Independent verification of judge claims before adoption: MMseqs2/CD-HIT
  binaries are not installable in this sandbox, so a documented 3-mer-Jaccard
  proxy with stated thresholds was locked instead; DRAMP/DBAASP public dumps
  confirmed freely downloadable (to be pinned with sha256 at fetch time).
