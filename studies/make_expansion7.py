"""Seventh expansion: extended base-model table, symbol glossary, study findings."""
import json, pathlib
R = pathlib.Path('results'); P = pathlib.Path('paper')

# app_basefull.tex: study11 all metrics for all models
v = json.load(open(R/'study11_veltri_stack.json'))
rows = []
names = {'base_cnn7':'PepCNN seed 7','base_cnn27':'PepCNN seed 27','base_cnn77':'PepCNN seed 77',
         'base_gnn11':'PepGNN seed 11','base_gnn111':'PepGNN seed 111','base_rf':'RF (descriptors)',
         'ensemble_avg_tuned':'Average ensemble','ensemble_stack_tuned':'Stacked ensemble'}
for k, label in names.items():
    m = v[k]
    thr = f", thr {m['threshold']:.3g}" if 'threshold' in m else ""
    rows.append(f"{label} & {100*m['sens']:.2f} & {100*m['spec']:.2f} & {100*m['acc']:.2f} & "
                f"{m['mcc']:.4f} & {100*m['auc']:.2f} & {100*m['auprc']:.2f} \\\\")
(P/'app_basefull.tex').write_text(r"""
\section{Complete study-11 metrics}\label{app:basefull}
Table~\ref{tab:basefull} gives every metric for every model of the Veltri
2018 study, including AUPRC, which the main text omits for space. The two
ensemble rows used the tune partition for their thresholds (0.56 average,
0.575 stack); base models use the default 0.5.

\begin{table}[h]\centering\small
\caption{All study-11 models on the Veltri test partition (percentages
except MCC).}
\label{tab:basefull}
\begin{tabular}{lcccccc}
\hline
Model & SENS & SPEC & ACC & MCC & AUC & AUPRC \\
\hline
""" + "\n".join(rows) + r"""
\hline
\end{tabular}
\end{table}

The AUPRC column reinforces the ranking: the stack (97.03) leads every
base model, and the GNN seeds trail on every metric, consistent with the
capacity and pooling argument of Appendix~\ref{app:arch}. Seed variance
among the CNNs is small on AUC (95.67--96.04) but larger on MCC
(0.759--0.805), the expected pattern when only the threshold varies.
""")

# app_symbols.tex: glossary
(P/'app_symbols.tex').write_text(r"""
\section{Symbol glossary}\label{app:symbols}
\begin{table}[h]\centering\small
\caption{Notation used throughout.}
\label{tab:symbols}
\begin{tabular}{ll}
\hline
Symbol & Meaning \\
\hline
$s = s_1\cdots s_n$ & peptide sequence, residues from the 20-letter alphabet \\
$L$ & padded sequence length (60 pilot, 150 large, 200 study 14) \\
$K_k(s)$ & set of contiguous $k$-mers of $s$ \\
$J_k$, $\nu_k$ & $k$-mer Jaccard similarity and novelty (Eq.~\eqref{eq:novelty}) \\
$Q(s;p)$ & expected net charge at pH $p$ (Eq.~\eqref{eq:netcharge}) \\
$\mu_H$, $\mu_i(s;w,\delta)$ & Eisenberg hydrophobic moment (Eq.~\eqref{eq:hmoment}) \\
$A$, $A_{\mathrm{norm}}$ & chain-window adjacency and its row normalization \\
$\hat p(s)$ & model posterior $P(\mathrm{AMP}\mid s)$ \\
$t^{*}$ & Bayes decision threshold (Eq.~\eqref{eq:bayesthr}) \\
$A$ (AUROC context) & Mann--Whitney AUROC statistic (Prop.~\ref{prop:mw}) \\
$Q_1, Q_2$ & placement probabilities in the Hanley--McNeil variance \\
$B$ & bootstrap or permutation resample count \\
$\alpha(n)$ & inverse Ackermann factor in union-find complexity \\
$\omega$ & positive-class weight in the weighted BCE \\
$c_{01}, c_{10}$ & false-positive and false-negative costs \\
$\mathrm{BS}$ & Brier score (Eq.~\eqref{eq:brier}) \\
\hline
\end{tabular}
\end{table}
""")

# extend app_studies with key findings narrative
(P/'app_findings.tex').write_text(r"""
\section{Key findings by study}\label{app:findings}
\begin{description}
\item[Study 1 (pilot).] On 2{,}000-reviewed-positive pilot data, PepCNN
reached AUROC 0.95 within minutes, establishing the pipeline end to end;
the annealing designer produced its first 24 candidates here.
\item[Study 7 (scale-up).] Cluster-aware splitting cut apparent
performance relative to a random split, the first signal that homology
inflates random-split estimates; the frozen screen ensemble dates from
this study and is never retrained afterward.
\item[Study 8 (screen).] 852 reviewed uncharacterized proteins scored;
41 pass the novelty floor; the identical P0ACW4/P0ACW5 sequence pair tops
the list.
\item[Study 9 (APD6).] The frozen ensemble scores APD6-2024 natural AMPs
high on average, a sanity check that the model prefers known AMPs to
uncharacterized proteins, without constituting validation of any
individual candidate.
\item[Study 10 (Veltri base).] On the authors' exact splits, PepCNN
reaches AUC 95.77, the RF 93.32, the GNN 88.18; composition carries much
of the benchmark signal.
\item[Study 11 (stack).] Stacking with tune-partition weights reaches ACC
91.22, MCC 0.8246, AUC 96.61: near-parity with the published AMP Scanner
row, behind ACEP; reported as descriptive, without a paired test.
\item[Study 12 (candidates).] Nine named candidates with novelty and
physicochemical audits and explicit falsifiable MIC predictions.
\item[Study 13 (filtered CV).] A Feb2020 CV that silently filtered records
produced numbers not comparable to the published ones; retained as an
ablation and superseded by study 14.
\item[Study 14 (full CV).] All 4{,}042 Feb2020 records, 10-fold CV: ACC
90.33, MCC 0.8076, AUC 96.36 versus published 89.9/0.799/96.2 ---
descriptive parity with a slight edge, no significance claim possible.
\item[Study 15 (ablation).] Removing the P0ACW4 signal peptide drops the
ensemble score from 0.9082 to 0.6299: a negative control showing the
score partly tracks trafficking signal.
\item[Study 16 (homology).] 15.98\% of held-out Feb2020 sequences share
an exact decamer with training; random-CV results are interpolation.
\item[Study 17 (strata).] The large-dataset test AUROC is 0.9206 overall,
0.9643 on the tiny reviewed-positive subgroup, 0.9121 on the confounded
majority: the headline number cannot be attributed to AMP function.
\item[Study 18 (composition).] First-order composition barely separates
the label $\times$ review-status groups; the learnable confound is mostly
a length artifact, directing the fix toward length-matched reviewed-only
redesign.
\end{description}
""")
print("expansion7 done")
