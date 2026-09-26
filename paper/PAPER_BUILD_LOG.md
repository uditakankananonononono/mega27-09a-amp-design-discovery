# Paper-build evidence log, 2026-09-27

Branch: `paper-build`, separate from the science builder's `main`. Paper paths only.

- Added locked follow-up for study 19 (`results/study19_feb2020_v2.json`): tenfold mean MCC 0.8328, accuracy 0.9161, AUROC 0.9688; recorded gates true. Published reference point estimates are not paired predictions on these exact folds.
- Added study 20 cluster-held-out stress test (`results/study20_cluster_holdout.json`): 4,001 clusters, MCC 0.8369 and AUROC 0.9701; sequence-cluster definition does not exclude remote homologs.
- Added study 21 ten-candidate triage (`results/study21_rescreen.json`): all named T1 under the specified reference audit, seven signal-peptide flags and three hemolysis-risk flags. No biological activity or safety was measured.
- Built `paper/main.pdf` with pdfLaTeX twice: 52 pages by `pdfinfo`. `mathptmx` maps body Times-style text to Nimbus Roman in this environment, **not licensed Times New Roman**; `pdffonts` confirms Nimbus Roman. Genuine TNR remains a font gate pending access to licensed fonts. The original 50-page PDF also used Nimbus Roman. Visually inspected final table and references; pending-work bullets corrected after the first preview exposed clipping.
