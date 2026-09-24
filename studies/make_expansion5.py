"""Fifth expansion: compositional audit section (study 18)."""
import json, pathlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
R = pathlib.Path('results'); P = pathlib.Path('paper')

d = json.load(open(R/'study18_compositional_audit.json'))
order = ['pos_reviewed', 'pos_unreviewed', 'neg_reviewed']
labels = ['pos\nreviewed', 'pos\nunreviewed', 'neg\nreviewed']
fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.0))
for ax, key, title in zip(axes, ['length_mean_sd', 'kd_mean_sd', 'net_charge_mean_sd'],
                          ['length (aa)', 'mean KD hydropathy', 'net charge (pH 7 proxy)']):
    m = [d[k][key][0] for k in order]; s = [d[k][key][1] for k in order]
    ax.bar(labels, m, yerr=s, capsize=3, color=['steelblue','lightsteelblue','indianred'])
    ax.set_title(title)
fig.tight_layout(); fig.savefig(P/'fig_composition.png', dpi=150); plt.close(fig)

rows = []
for k in order:
    v = d[k]
    rows.append(f"{k.replace('_',' ')} & {v['n']:,} & {v['length_mean_sd'][0]:.1f}$\\pm${v['length_mean_sd'][1]:.1f} & "
                f"{v['kd_mean_sd'][0]:+.2f}$\\pm${v['kd_mean_sd'][1]:.2f} & "
                f"{v['net_charge_mean_sd'][0]:+.2f}$\\pm${v['net_charge_mean_sd'][1]:.2f} & "
                f"{100*v['aromatic_frac_mean_sd'][0]:.1f} \\\\")
(P/'compositional_audit.tex').write_text(r"""
\subsection{Compositional audit of the confound}\label{sec:compaudit}
What separates the classes if not labels? Study 18 computes first-order
composition for the three label $\times$ review-status groups of the large
dataset (Table~\ref{tab:compaudit}, Figure~\ref{fig:composition}).

\begin{table}[h]\centering\small
\caption{Composition by label and review status over all 49{,}256 large
dataset records (train+val+test). Charge is the integer formal-charge
proxy of Appendix~\ref{app:encoding}.}
\label{tab:compaudit}
\begin{tabular}{lrcccc}
\hline
Group & $n$ & Length & KD hydropathy & Net charge & Aromatic \% \\
\hline
""" + "\n".join(rows) + r"""
\hline
\end{tabular}
\end{table}

\begin{figure}[h]\centering
\includegraphics[width=\textwidth]{fig_composition.png}
\caption{Length, hydropathy and charge by group (bars: mean $\pm$ sd).}
\label{fig:composition}
\end{figure}

The result is a caution against a simple story. First-order composition
barely separates the groups: mean charge spans only $+2.4$ to $+3.2$, mean
hydropathy $-0.29$ to $+0.03$. The one large gap is length: reviewed
positives average 66 aa against 93--97 aa for both other groups, a
curation artifact (reviewed AMP records skew to mature short peptides,
unreviewed records to longer precursors). Length is trivially learnable by
a CNN through padding statistics, so part of the 0.92 AUROC plausibly
reflects a length prior rather than antimicrobial chemistry; the rest must
come from subtler annotation-correlated sequence patterns (source
organisms, precursor context, pipeline-specific motifs). This is exactly
why Section~\ref{sec:labelconfound} retires the large-dataset headline:
the model learned something real, but what it learned cannot be cleanly
attributed to AMP function, and a length-matched reviewed-only redesign is
the correct fix, recorded as future work.
""")
print("expansion5 done")
