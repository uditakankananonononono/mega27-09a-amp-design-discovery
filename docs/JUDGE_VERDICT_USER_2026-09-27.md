# USER-PROVIDED JUDGE VERDICT - 2026-09-27 (THE counted judge round for 09a)

Source: WhatsApp 2026-09-27 11:26:52 IST,
wamid.HBgMOTE4MTM0MDk4NTcxFQIAEhgWM0VCMEQ0OUM1NzUyMTgxOTJDQzE2NwA=,
verified author=user (5,875 body bytes). Pasted critique (ChatGPT-sourced per
the courier route) PROVIDED BY THE USER; under her 2026-09-27 rule
("EACH PROJECTS NEED ONE FROM ME TO PASS") this verdict is 09a's ONE counted
judge round. Content below is verbatim, unedited. Factual claims in it are
verified independently before adoption, as always.

---

20 Weaknesses
#	Weakness	Why It Matters for ISEF CBIO
1	Label-provenance confound invalidates the main scale-up result	Positives are 88% unreviewed (TrEMBL) while all negatives are reviewed (Swiss-Prot). The model can learn review status, not AMP function. AUROC 0.921 is not a clean recognition measure.
2	External APD6 validation still uses reviewed negatives	The confound persists: APD6 positives are validated AMPs, but negatives are the same reviewed UniProt set. The AUROC 0.926 cannot be interpreted as cross-database generalization.
3	The GNN underperforms the CNN at every scale	Chain-window adjacency is a weak prior; message passing adds parameters without adding information. This is reported as a negative architectural result, but it weakens the “two architectures” framing.
4	No wet-lab validation of any candidate	All discovery candidates are in silico. ISEF CBIO allows computational projects, but the lack of any experimental follow-up limits impact.
5	Signal-peptide ablation shows the model may learn trafficking signals	Removing the signal peptide drops the ensemble score from 0.908 to 0.630. The high score partly reflects secretion, not antimicrobial chemistry.
6	Published benchmark comparisons are not statistically significant	Near-parity with AMP Scanner, below ACEP. No paired significance test possible. Descriptive parity is a weak claim.
7	Random-split homology leakage in published benchmarks	15.98% of held-out Feb2020 sequences share an exact decamer with training. Cross-validation estimates interpolation, not extrapolation.
8	Negative class is “not annotated antimicrobial”	Not the same as experimentally verified non-AMP. Some negatives may be unrecognized AMPs.
9	Pilot study uses synthetic identifiers	Provenance unverifiable; pilot is a software sanity check, not a benchmark.
10	10-mer clustering is a weak homology control	Two sequences sharing no exact decamer can still be homologs. Alignment-based clustering (MMseqs2) is needed.
11	Annealing designer optimizes the same model that scores it	Echo chamber: generated sequences are scored by the ensemble that guided generation. No independent validation.
12	Novelty metric is 3-mer Jaccard	Weak; high novelty does not imply absence from databases or biological novelty.
13	No activity prediction, only binary classification	Cannot predict MIC or potency. Discovery candidates are ranked hypotheses, not calibrated probabilities.
14	Compute budget limits model size and epochs	2 CPU, 2 GB RAM; 30k training subsample, 8 epochs. Larger models might change results.
15	The paper is 56 pages, 34 tables, 18 appendices	ISEF judges have 12 minutes. The signal-to-noise ratio is poor.
16	The main positive result is confounded	The CNN AUROC 0.921 and APD6 AUROC 0.926 are measured on confounded data. The clean reviewed-only subgroup (239 positives) is too small for robust claims.
17	No family-level analysis	The model may miss entire AMP classes (e.g., anionic, disulfide-rich). Sensitivity 0.601 at 5% FPR suggests many false negatives.
18	The ensemble is a simple average	No learned weights or calibration. The stack on Veltri uses logistic regression, but the discovery screen uses unweighted mean.
19	The “discovery” screen is over 852 uncharacterized proteins	Top hits are cysteine-rich, hydrophobic, signal-peptide-like. This may reflect the confounded training set, not true AMPs.
20	No pre-registration of the discovery screen	Candidates were selected after looking at model scores. Selection is conditional; scores are not unbiased validation metrics.
20 Additions / Improvements
#	Addition	How It Strengthens the Project
1	Fetch reviewed-only positives and match negatives on review status and organism	Removes the confound. Rebuild splits, retrain, re-evaluate. This is the top priority.
2	Use DBAASP activity data for potency regression	Replaces binary labels with continuous MIC. Makes predictions clinically meaningful.
3	Replace 10-mer clustering with MMseqs2 alignment-based clustering	Stricter homology control; better generalization estimates.
4	Add structure-aware GNN using predicted contact maps or ESM embeddings	Tests whether graph models help when the graph carries real information.
5	Collaborate with a wet lab to test top 3–5 candidates	Even a small MIC assay would transform the project from in silico to validated.
6	Pre-register the discovery screen before scoring	Locks the hypothesis and prevents post-hoc selection.
7	Report per-family performance	Show which AMP classes are missed (e.g., anionic, disulfide-rich).
8	Add calibration (Platt scaling, isotonic regression)	Makes ensemble outputs interpretable as probabilities.
9	Optimize ensemble weights on a validation set	Better than simple average.
10	Use a held-out external cohort not from UniProt/APD6	E.g., DBAASP or a recently published AMP set.
11	Add negative controls: scrambled sequences, known non-AMPs	Tests whether model scores are specific.
12	Ablate encoding features	Which channels (one-hot vs physicochemical) drive performance?
13	Test length confound explicitly	Train on length-matched reviewed-only data.
14	Use a stricter novelty metric (alignment-based)	E.g., BLAST or HMMER against all public databases.
15	Provide a simple web tool for screening	Increases usability and impact.
16	Shorten the paper to 12 slides for ISEF	Focus on: confound discovery, CNN performance, APD6 validation, top candidates, wet-lab protocol.
17	Lead with the label-provenance audit as the primary contribution	It is a methodological warning for the entire AMP prediction field.
18	Compare to more recent AMP predictors (AMPlify, Macrel, etc.)	Positions the work against current state of the art.
19	Propose a benchmark standard	Report review-status cross-tabs, length matching, cluster splits.
20	Add a decision tree for practitioners	“If you use UniProt KW-0929, here is how to avoid the confound.”
