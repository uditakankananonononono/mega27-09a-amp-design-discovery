"""Second expansion: figures + dataset glossary + study walkthrough + CLI appendix."""
import json, pathlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R = pathlib.Path('results'); P = pathlib.Path('paper')

# Fig: study14 per-fold metrics vs published means
d = json.load(open(R/'study14_feb2020_full.json'))
folds = [f['fold'] for f in d['folds']]
fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.2))
for ax, key, pub in zip(axes, ['acc','mcc','auc'], [0.899, 0.799, 0.962]):
    vals = [f[key] for f in d['folds']]
    ax.plot(folds, vals, 'o-', label='this work (per fold)')
    ax.axhline(pub, color='r', ls='--', label=f'published mean {pub}')
    m, s = d['mean_sd'][key]
    ax.axhline(m, color='k', ls=':', label=f'our mean {m:.3f}')
    ax.set_xlabel('fold'); ax.set_title(key.upper()); ax.set_ylim(min(min(vals),pub)-0.02, max(max(vals),pub)+0.02)
axes[0].legend(fontsize=7, loc='lower right')
fig.tight_layout(); fig.savefig(P/'fig_feb2020_folds.png', dpi=150); plt.close(fig)

# Fig: per-fold homology overlap
h = json.load(open(R/'study16_homology_audit.json'))
fig, ax = plt.subplots(figsize=(6.5, 3.0))
ax.bar([f['fold'] for f in h['folds']], [100*f['fraction'] for f in h['folds']], color='steelblue')
ax.axhline(100*h['overall_fraction'], color='r', ls='--', label=f"overall {100*h['overall_fraction']:.1f}%")
ax.set_xlabel('fold'); ax.set_ylabel('% held-out sharing an exact decamer'); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(P/'fig_homology_folds.png', dpi=150); plt.close(fig)

# app_datasets.tex
(P/'app_datasets.tex').write_text(r"""
\section{Dataset glossary}\label{app:datasets}
\begin{table}[h]\centering\scriptsize
\caption{Every dataset used in this work, its role, and its caveats.}
\label{tab:datasets}
\begin{tabular}{p{3.2cm}p{1.6cm}p{2.6cm}p{3.4cm}p{4.2cm}}
\hline
Dataset & Size & Role & Source & Caveat \\
\hline
UniProt KW-0929 positives (10--150 aa) & 27{,}392 fetched; 24{,}628 used &
Positive class, large benchmark & UniProt REST, keyword KW-0929 &
21{,}742 of 24{,}628 used records are unreviewed (Section~\ref{sec:labelconfound}) \\
UniProt reviewed non-AMP pool & 105{,}041 fetched; 24{,}628 used &
Negative class, length-matched & UniProt REST, reviewed, NOT KW-0929 &
All reviewed: review-status asymmetry vs positives \\
APD6 natural AMPs (2024) & 3{,}306 & External validation corpus &
aps.unmc.edu downloads & Post-2024 list; successor of APD3 \\
UniProt uncharacterized (10--60 aa) & 852 & Discovery screen pool &
UniProt REST, reviewed, uncharacterized & Scores are hypotheses, not activities \\
Veltri 2018 original splits & 1{,}424 / 708 / 1{,}424 &
Published benchmark (train/tune/test) & authors' GitHub &
Fixed single split; identity-reduced by authors \\
Feb2020 production dataset & 2{,}021 + 2{,}021 &
Published benchmark (10-fold CV) & authors' site &
Published fold assignments unavailable; 15.98\% decamer overlap
(Appendix~\ref{app:homology}) \\
\hline
\end{tabular}
\end{table}
""")

# app_studies.tex: walkthrough
(P/'app_studies.tex').write_text(r"""
\section{Study-by-study reproducibility walkthrough}\label{app:studies}
Every study is a single script under \texttt{studies/} writing one JSON
(and figures where noted) under \texttt{results/}. Commands are run from
the repository root; hermetic fixtures cover network access in tests.

\begin{table}[h]\centering\scriptsize
\caption{Study inventory. Runtimes are sandbox wall-clock where measured.}
\label{tab:studies}
\begin{tabular}{p{1.3cm}p{6.2cm}p{4.6cm}p{1.8cm}}
\hline
Script & Purpose & Output & Runtime \\
\hline
\texttt{study01} & Pilot benchmark: reviewed KW-0929 (8--60 aa) vs
length-matched reviewed non-AMP; PepCNN/PepGNN/RF/logreg; annealing
designer & \texttt{study01\_amp.json}, ROC figure & minutes \\
\texttt{study07} & Scale-up on homology-aware clustered split of the
24{,}628 + 24{,}628 large dataset; trains the frozen screen ensemble &
\texttt{study07\_scaleup.json}, checkpoints & hours \\
\texttt{study08} & Discovery screen: score 852 reviewed uncharacterized
proteins with the frozen ensemble & \texttt{study08\_discovery.json} & minutes \\
\texttt{study09} & External validation: score APD6 2024 natural AMPs with
the same frozen ensemble & \texttt{study09\_apd\_validation.json} & minutes \\
\texttt{study10} & Base-model head-to-head on the Veltri 2018 splits &
\texttt{study10\_veltri.json} & $\sim$1 h \\
\texttt{study11} & Stacked ensemble on the Veltri 2018 splits, tuned on
the tune partition only & \texttt{study11\_veltri\_stack.json} & $\sim$2 h \\
\texttt{study12} & Candidate naming, novelty audit (APD6, KW-0929),
physicochemical audit, falsifiable predictions &
\texttt{study12\_candidates.json} & minutes \\
\texttt{study13} & Filtered-record Feb2020 CV (retained as ablation) &
\texttt{study13\_feb2020\_cv.json} & $\sim$15 min \\
\texttt{study14} & Full-record Feb2020 10-fold CV, all 4{,}042 sequences &
\texttt{study14\_feb2020\_full.json} & $\sim$15 min \\
\texttt{study15} & P0ACW4 precursor/signal/mature segment ablation &
\texttt{study15\_signal\_ablation.json} & seconds \\
\texttt{audit\_homology} & Exact-decamer train/test overlap per fold &
\texttt{study16\_homology\_audit.json} & seconds \\
\texttt{study17} & Post-hoc review-status stratification of the large
test set & \texttt{study17\_review\_strata.json} & minutes \\
\hline
\end{tabular}
\end{table}

The full pipeline is: \texttt{fetch\_large.py} and
\texttt{build\_large\_dataset.py} (data, with provenance JSON),
\texttt{fetch\_published\_benchmarks.py} (pinned external benchmarks,
Appendix~\ref{app:sources}), the study scripts above in numerical order,
then the paper-asset generators (\texttt{make\_paper\_assets.py},
\texttt{make\_sota\_section.py}, \texttt{make\_figs2.py},
\texttt{make\_appendices2.py}, \texttt{make\_expansion2.py}).
""")

# app_cli.tex
(P/'app_cli.tex').write_text(r"""
\section{Screening CLI walkthrough}\label{app:cli}
The repository ships a screening command, \texttt{pepdesign screen},
wrapping the frozen study-07 ensemble. A session:

\begin{verbatim}
$ pepdesign screen --input uncharacterized.fasta --top 25 \
    --novelty-floor 0.90 --out screen.json
screened 852 sequences; 41 pass the novelty floor
top hit: sp|P0ACW4|YDCA_ECOLI score 0.9082 novelty 0.938
\end{verbatim}

The JSON output records, per sequence, the ensemble score, both sub-scores,
the 3-mer novelty against the reference corpus, and whether the novelty
floor passed. The bundled demo (\texttt{results/tool\_screen\_demo.json})
screens the same 852-protein pool as study 08 and reproduces its top hits,
including the duplicate-accession pair P0ACW4/P0ACW5 (identical 57-residue
sequence in two species records), which the tool reports verbatim rather
than silently merging.

Design use: the annealing designer is exposed as
\texttt{pepdesign design --seeds 24 --steps 250}, reproducing the study-07
anneal batch; the scoring function, cooling schedule and novelty floor are
the ones derived in Section~\ref{sec:results} and Appendix~\ref{app:hyper}.
""")
print("expansion2 written")
