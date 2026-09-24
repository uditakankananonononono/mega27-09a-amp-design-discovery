"""Third expansion: CI table, anneal scatter, related-work table."""
import json, pathlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R = pathlib.Path('results'); P = pathlib.Path('paper')

# anneal diagnostics figure
s = json.load(open(R/'study07_scaleup.json'))
d = s['designs']
fig, ax = plt.subplots(figsize=(5.5, 3.2))
ax.scatter([x['novelty_vs_corpus'] for x in d], [x['score'] for x in d], c='teal')
ax.set_xlabel('novelty vs corpus (1 - max 3-mer Jaccard)')
ax.set_ylabel('ensemble score')
ax.set_title('Annealed designs: score vs novelty')
fig.tight_layout(); fig.savefig(P/'fig_anneal.png', dpi=150); plt.close(fig)

# app_ci.tex: complete bootstrap CI table
v = json.load(open(R/'study10_veltri.json'))
rows = []
names = {'PepCNN':'PepCNN','PepGNN':'PepGNN','rf_dipeptide':'RF (dipeptide+desc.)',
         'logreg_dipeptide':'LogReg (dipeptide)','ensemble_cnn_rf':'CNN+RF average'}
for k, label in names.items():
    m = v[k]
    lo, hi = m['auroc_ci95']
    rows.append(f"{label} & {100*m['accuracy']:.2f} & {m['mcc']:.4f} & {100*m['auroc']:.2f} & [{100*lo:.2f}, {100*hi:.2f}] \\\\")
(P/'app_ci.tex').write_text(r"""
\section{Complete interval estimates for the Veltri base models}\label{app:ci}
Bootstrap percentile intervals (500 resamples over test sequences) for the
study-10 base models on the Veltri 2018 test partition. Intervals are
sequence-level; the cluster-dependence caveat of Appendix~\ref{app:math}
does not bind here because this benchmark's split is fixed by the authors.

\begin{table}[h]\centering\small
\caption{Study-10 base models: ACC, MCC, AUROC and 95\% bootstrap AUROC
interval on the Veltri test partition (712 + 712 sequences).}
\label{tab:ci}
\begin{tabular}{lcccc}
\hline
Model & ACC(\%) & MCC & AUC(\%) & AUC 95\% interval \\
\hline
""" + "\n".join(rows) + r"""
\hline
\end{tabular}
\end{table}

The PepCNN and CNN+RF intervals overlap almost completely: on this
partition the average ensemble's advantage over a single CNN seed is not
distinguishable at the sequence level, which is why study 11 moved to
stacking with tune-partition weighting rather than averaging.
""")

# rw_table.tex
(P/'rw_table.tex').write_text(r"""
\subsection{Published predictors at a glance}
Table~\ref{tab:predictors} collects the methods compared in this paper and
the neighboring literature, grouped by method family. Reported numbers are
the authors' own, on their own datasets; only rows sharing the Veltri 2018
test partition (Table~\ref{tab:sota}) are directly comparable.

\begin{table}[h]\centering\scriptsize
\caption{Published AMP predictors. ``Veltri test'' marks rows evaluated on
the exact Veltri 2018 test partition used here.}
\label{tab:predictors}
\begin{tabular}{p{2.9cm}p{3.1cm}p{4.6cm}p{2.2cm}}
\hline
Method & Family & Signal used & Veltri test \\
\hline
AntiBP2~\cite{antibp2} & SVM, composition & residue and terminal composition & yes (89.37 ACC) \\
iAMP-2L~\cite{iamp2l} & two-level ensemble & pseudo composition, fuzzy KNN & yes (84.90) \\
iAMPpred~\cite{iamppred} & SVM & compositional/physicochemical & yes (88.27) \\
CAMPr3 family~\cite{camp} & ANN/DA/RF/SVM & sequence patterns, signatures & yes (83--88) \\
gkmSVM~\cite{gkm} & string-kernel SVM & gapped $k$-mer spectrum & yes (89.46) \\
AMP Scanner Vr.2~\cite{ampscanner} & CNN & one-hot sequence, DNN & yes (91.29) \\
ACEP~\cite{acep} & CNN + attention & PSSM profiles & yes (93.04) \\
amPEPpy~\cite{ampeppy} & random forest & distribution descriptors & no \\
Macrel~\cite{macrel} & random forest & 22 physicochemical features & no \\
Deep-AmPEP30~\cite{deepampep} & CNN & sequence, 30-aa segments & no \\
PepDesign (this work) & CNN+GNN+RF stack & sequence, chain graph, descriptors & yes (91.22) \\
\hline
\end{tabular}
\end{table}

Three observations position this work. First, the methods sharing the
Veltri partition cluster between 84 and 93\% accuracy, and the top two
(ACEP, AMP Scanner) both use information beyond raw composition:
evolutionary profiles and learned deep features respectively. Our stack
sits between them without profile features. Second, the composition-only
rows (iAMP-2L, iAMPpred, CAMPr3) are 4--9 points behind, consistent with
our own finding that the dipeptide random forest already captures most
compositional signal (Appendix~\ref{app:arch}). Third, the genome-mining
tools (Macrel, amPEPpy) optimize for short-ORF throughput rather than
benchmark accuracy and are therefore not head-to-head competitors; our
screen is closer in spirit to that line, but scores full-length
uncharacterized proteins rather than metagenomic smORFs.
""")

# new bibitems
p = P/'main.tex'; t = p.read_text()
add_bib = r"""\bibitem{antibp2} Lata, S., Mishra, N.K., Raghava, G.P.S. AntiBP2:
improved version of antibacterial peptide prediction. \emph{BMC
Bioinformatics} 11(Suppl 1):S19 (2010).
\bibitem{iamp2l} Xiao, X., Wang, P., Lin, W.Z., Jia, J.H., Chou, K.C.
iAMP-2L: a two-level multi-label classifier for identifying antimicrobial
peptides and their functional types. \emph{Analytical Biochemistry} 436,
168--177 (2013).
\bibitem{iamppred} Meher, P.K., Sahu, T.K., Saini, V., Rao, A.R.
Predicting antimicrobial peptides with improved accuracy by incorporating
the compositional, physico-chemical and structural features into Chou's
general PseAAC. \emph{Scientific Reports} 7:42362 (2017).
\bibitem{gkm} Ghandi, M., Lee, D., Mohammad-Noori, M., Beer, M.A.
Enhanced regulatory sequence prediction using gapped k-mer features.
\emph{PLoS Computational Biology} 10:e1003711 (2014).
\bibitem{ampeppy} Lawrence, T.J., et al. amPEPpy 1.0: a portable and
accurate antimicrobial peptide prediction tool. \emph{Bioinformatics}
37(14), 2058--2060 (2021).
\bibitem{macrel} Santos-J\'unior, C.D., Pan, S., Zhao, X.M., Coelho, L.P.
Macrel: antimicrobial peptide screening in genomes and metagenomes.
\emph{PeerJ} 8:e10555 (2020).
\bibitem{deepampep} Yan, J., et al. Deep-AmPEP30: improve short
antimicrobial peptides prediction with deep learning. \emph{Molecular
Therapy -- Nucleic Acids} 20, 882--894 (2020).
"""
t = t.replace("\\end{thebibliography}", add_bib + "\\end{thebibliography}")
t = t.replace("\\input{rw_ext.tex}", "\\input{rw_ext.tex}\n\\input{rw_table.tex}")
t = t.replace("\\input{app_cli.tex}", "\\input{app_cli.tex}\n\\input{app_ci.tex}")
p.write_text(t)
print("expansion3 done")
