# PREREGISTRATION ADDENDUM 6 - Q3 DBAASP potency regression (locked 2026-09-27, before any label pull or training)

Source: addendum-4 Q3 (=A2), itself from the 2026-09-26 master spec. This addendum
locks the dataset construction, target transform, model, protocol, and metrics
BEFORE any MIC labels are pulled or any model is trained.

## Data gap disclosure (resolved by this addendum)
The sha256-verified in-repo DBAASP copy (data/raw/novelty_refs/dbaasp_all.fasta,
manifest sha256 fbaebb52...) contains SEQUENCES ONLY - no activity values. MIC
labels are therefore pulled live from the public DBAASP API
(https://dbaasp.org/peptides/<id>), one record per in-repo canonical sequence,
matched by the DBAASPR_<id> header. Raw API responses are stored under
data/raw/dbaasp_api/ (gitignored) with a sha256 manifest committed to
docs/DBAASP_API_MANIFEST_2026-09-27.md. No label is read or inspected before
this lock; the pull script writes labels to disk without summarizing them.

## Locked dataset construction
- Universe: the 546 canonical 5-150 aa monomer sequences already in-repo
  (dbaasp_all.fasta), deduplicated by exact sequence string.
- Assay filter: targetActivities entries with activityMeasureGroup.name == "MIC"
  ONLY. Hemolytic/antibiofilm/other activity types excluded.
- Primary cohort: targetSpecies.name == "Escherichia coli ATCC 25922", unit
  µg/ml only (records in µM are converted ONLY when the peptide molecular weight
  is derivable from the plain 20-aa sequence; conversion rule
  µg/ml = µM x MW(Da) / 1000, MW via average residue masses; peptides with
  non-standard residues stay µg/ml-only or are dropped from the primary cohort).
  Excluded-conversion count is reported.
- Replicates: median of all qualifying MIC values per peptide.
- Target: y = log10(median MIC, µg/ml). Locked.
- Secondary (exploratory only, clearly labeled, no headline): all-species pooled
  µg/ml cohort, same rules. No success claims from secondary cohorts.

## Locked model and protocol
- Model: HistGradientBoostingRegressor on the same physicochemical feature set
  used by the HGB/RF components of the study19 ensemble. The GNN component is
  classification-only and is EXCLUDED (disclosed, not a silent omission).
- Hyperparameters: library defaults except max_iter=300, random_state=11. Locked.
- Splits: 10-fold cluster-held-out CV, MMseqs2 clusters recomputed on THIS
  sequence set with the locked Q2 params (min-seq-id 0.3, cov 0.8, cov-mode 0,
  union-find single-linkage), rng(11) greedy fold assignment, same as study23.
- Null model: global-median predictor on the same folds.
- Falsifier: one label-shuffle run (rng 11); expectation rho ~= 0. A falsifier
  rho > 0.2 invalidates the arm and is reported as a failed control.

## Locked metrics (primary first)
1. Spearman rho (pooled out-of-fold predictions) - primary.
2. MAE in log10(µg/ml) units.
3. R^2.
All reported with the null-model deltas. ANY direction of result is reported
as-is; there is no success threshold and no claim that potency prediction beats
a published benchmark unless the numbers say so on their face.

## Deviations
None at lock time. Any deviation discovered during execution is disclosed
in-doc here and in the results JSON, never silently.

## DEVIATION-1 (disclosed 2026-09-27, before any label fetch)
1. The in-repo dbaasp_all.fasta is DEGENERATE: 546 FASTA entries but only 12
   unique sequences (original novelty-ref pull bug; headers carry parent
   dbaaspId repeated per monomer chain). Header-ID matching is impossible.
2. Consequence for the locked universe: the "546 canonical sequences already
   in-repo" do not exist as distinct sequences. Universe is redefined as ALL
   canonical 20-aa 5-150 aa monomer sequences from the live DBAASP catalog
   (25,542 records, paginated list endpoint carries sequences only, no
   activity data). The label pull (full /peptides/<id> records) begins only
   after this disclosure; the blind-pull property is preserved.
3. A corrected novelty reference (data/raw/novelty_refs/
   dbaasp_catalog_2026-09-27.fasta, sha256 in the pull manifest) replaces the
   degenerate file. Impact on study21's DBAASP novelty arm (which read the
   degenerate file) is being re-audited and will be corrected in the open.
4. Everything else in this addendum (assay filter, primary cohort, target
   transform, model, folds, metrics, falsifier) is unchanged.

## CORRECTION-1 to DEVIATION-1 (disclosed 2026-09-27, ~14:51 IST)
The DEVIATION-1 timeline statement "the pull script writes labels to disk
without summarizing them / the label pull begins only after this disclosure"
was imprecise. Exact timeline:
1. ~14:14:58-14:15:01, BEFORE addendum 6 was locked: two manual API-structure
   probes (curl /peptides/1 and /peptides/3) printed one MIC record each to
   the terminal (peptide 1: MIC 14.5 uM, E. coli ATCC 25922; peptide 3:
   MIC 16 ug/ml, S. aureus MRSA). These two label values were visible before
   the lock. They played no role in any locked choice - target transform,
   model, folds, and metrics were inherited unchanged from addendum-4 Q3,
   which predates all contact with the API.
2. ~14:15:22, AFTER addendum 6 (59916c1) but BEFORE DEVIATION-1 (87b4c7f):
   the first scripted pull wrote 10 full records to disk (ids 1,3,4,6,9,18,
   19,21,22 + 1 more), unread and unsummarized (log shows counts only).
   Those files remain as raw data and are a subset of the redefined universe.
No label-informed decision occurred at any point; all other labels remain
uninspected until the locked pipeline runs. This correction is also being
reported to the parent for the record.
