# Paper-build evidence log, 2026-09-27

Branch: `paper-build`, separate from the science builder's `main`. Paper paths only.

- Added locked follow-up for study 19 (`results/study19_feb2020_v2.json`): tenfold mean MCC 0.8328, accuracy 0.9161, AUROC 0.9688; recorded gates true. Published reference point estimates are not paired predictions on these exact folds.
- Added study 20 cluster-held-out stress test (`results/study20_cluster_holdout.json`): 4,001 clusters, MCC 0.8369 and AUROC 0.9701; sequence-cluster definition does not exclude remote homologs.
- Added study 21 ten-candidate triage (`results/study21_rescreen.json`): all named T1 under the specified reference audit, seven signal-peptide flags and three hemolysis-risk flags. No biological activity or safety was measured.
- Built `paper/main.pdf` with pdfLaTeX twice: 52 pages by `pdfinfo`. `mathptmx` maps body Times-style text to Nimbus Roman in this environment, **not licensed Times New Roman**; `pdffonts` confirms Nimbus Roman. Genuine TNR remains a font gate pending access to licensed fonts. The original 50-page PDF also used Nimbus Roman. Visually inspected final table and references; pending-work bullets corrected after the first preview exposed clipping.

## Preregistered t=0.40 primary arm, 2026-09-27

- Added the committed study21 harder-split primary result from main commit 8413a88, without merging or editing science code/results. The preregistration is `docs/PREREG_ADDENDUM3_2026-09-27.md`; diagnostic is `results/study21_clusters_t0.4.json` (3,721 clusters, 3,536 singletons, largest 12). Main metric artifact is `results/study21_hard_split_t0.4.json`: ten folds, mean MCC 0.8337, SD 0.0198, lower one-sided 95% 0.8223; both fixed point-value gates pass. The 0.0009 margin over the stronger 0.8328 gate is narrow and does not establish external family generalization. The secondary t=0.35 and t=0.45 arms are pending and omitted.
- pdfLaTeX twice yielded 54 pages (from 52), Nimbus Roman rather than genuine Times New Roman. Visually checked pages 52-54 and shortened a source path that initially overflowed the page margin.

## Judge requirement amended, 2026-09-27 10:00 IST

The owner said "NOT 10 ROUNDS OOF CHATGPT CHECK JUST ONE WHICH I PROVIDE OK?" (authenticated WhatsApp message `wamid.HBgMOTE4MTM0MDk4NTcxFQIAEhgWM0VCMDJCMTZGRTVEMkQwMTFBQzc4MQA=`, 10:00:07 IST). For this project, the paper branch therefore marks the counted judge gate **0 of 1, PENDING her personally provided verdict**. A round initiated by agents, even through her ChatGPT account, remains historical or supplementary and does not meet the gate. Historical ten-round language in the science ledgers is not erased by this note. A project-specific user-pasted verdict must be traced and evaluated before completion is recorded. Supplementary Gemini/LLM consults do not count. No scientific result, page count or font gate changes here.

## Preregistered t=0.35 sensitivity arm, 2026-09-27

- Added committed `results/study21_hard_split_t0.35.json` from main c296ce3 and diagnostics from main 5e7f660. Mean MCC 0.8231 passes the published 0.799 point-value hard gate and fails the stronger 0.8328 gate. The one-sided lower 95% bound is 0.7972, below 0.799. Fold SD is 0.0447; the 0.35 diagnostic has 3,624 clusters, and its 50-element leading-size list cannot establish a singleton count. The primary t=0.40 arm stays primary; t=0.45 remains pending.
- Rebuilt twice: 56 pages (up from 54) in Nimbus Roman, not genuine TNR. Visual check pages 53-55 verified the new two-page section's readable fold and summary tables, with no clipping; references continue on page 56. The second sensitivity arm does not validate the discovery candidates or change the pending user-provided judge-verdict gate.

## Authorship-attribution cleanup, 2026-09-27 11:14 IST

The owner requested removal of the assistant's attribution from the papers (WhatsApp `wamid.HBgMOTE4MTM0MDk4NTcxFQIAEhgWM0VCMEY5MzY4M0Q4OUYwNjg4ODZDNwA=`). Removed agent/program-style byline and credit text from the editable paper source and PDF display, without substituting an author. Udita's own byline in 09b was preserved, with only the Instinct pipeline parenthetical removed. Manuscript PDF author metadata is empty. Literature references to other studies' authors and technical uses of "author numbering" are not authorship credits for this paper.
