"""Generate paper/sota_benchmark.tex from studies 11 and 14 JSONs."""
import json, pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
s11 = json.load(open(ROOT/"results/study11_veltri_stack.json"))
s13 = json.load(open(ROOT/"results/study14_feb2020_full.json"))

def pct(x): return f"{100*x:.2f}"

published_rows = [
 ("AntiBP2", .8791, .9080, .8937, .7876, .8936),
 ("CAMPr3-ANN", .8300, .8511, .8405, .6813, .8405),
 ("CAMPr3-DA", .8707, .8075, .8391, .6797, .8997),
 ("CAMPr3-RF", .9269, .8244, .8757, .7553, .9363),
 ("CAMPr3-SVM", .8862, .8047, .8455, .6933, .9062),
 ("iAMP-2L", .8399, .8586, .8490, .6983, .8490),
 ("iAMPpred", .8933, .8722, .8827, .7656, .9444),
 ("gkmSVM", .8834, .9059, .8946, .7895, .9498),
 ("AMPScanner (Veltri 2018)", .8988, .9269, .9129, .8261, .9630),
 ("ACEP (Fu 2020)", .9241, .9367, .9304, .8610, .9778),
]
st = s11["ensemble_stack_tuned"]
published_rows.append(("PepDesign stack (this work)", st["sens"], st["spec"], st["acc"], st["mcc"], st["auc"]))
tab1 = "\n".join(f"{n} & {pct(a)} & {pct(b)} & {pct(c)} & {d:.4f} & {pct(e)} \\\\" for n,a,b,c,d,e in published_rows)

base_rows = []
for k, v in s11.items():
    if k.startswith("base_"):
        base_rows.append(f"\\texttt{{{k[5:]}}} & {pct(v['acc'])} & {v['mcc']:.4f} & {pct(v['auc'])} \\\\")
avg = s11["ensemble_avg_tuned"]
base_rows.append(f"average (thr {avg['threshold']:.2f}) & {pct(avg['acc'])} & {avg['mcc']:.4f} & {pct(avg['auc'])} \\\\")
base_rows.append(f"\\textbf{{stack (thr {st['threshold']:.2f})}} & \\textbf{{{pct(st['acc'])}}} & \\textbf{{{st['mcc']:.4f}}} & \\textbf{{{pct(st['auc'])}}} \\\\")
tab2 = "\n".join(base_rows)

fold_rows = "\n".join(
    f"{r['fold']} & {pct(r['sens'])} & {pct(r['spec'])} & {pct(r['acc'])} & {r['mcc']:.4f} & {pct(r['auc'])} \\\\"
    for r in s13["folds"])
ms = s13["mean_sd"]; pl = s13["pooled_oof"]
fold_rows += ("\n\\hline\n\\multicolumn{6}{l}{mean $\\pm$ sd: " +
  f"SENS {pct(ms['sens'][0])}$\\pm${100*ms['sens'][1]:.1f}, SPEC {pct(ms['spec'][0])}$\\pm${100*ms['spec'][1]:.1f}, " +
  f"ACC {pct(ms['acc'][0])}$\\pm${100*ms['acc'][1]:.1f}, MCC {ms['mcc'][0]:.4f}$\\pm${ms['mcc'][1]:.3f}, " +
  f"AUC {pct(ms['auc'][0])}$\\pm${100*ms['auc'][1]:.1f}" + "} \\\\")
fold_rows += (f"\npooled OOF & {pct(pl['sens'])} & {pct(pl['spec'])} & {pct(pl['acc'])} & {pl['mcc']:.4f} & {pct(pl['auc'])} \\\\")
v = s13["verdict"]
verdict_tex = "\n".join(
    f"{k.upper()} & {pct(x['published'])} & {pct(x['ours_mean'])} & {pct(x['ours_pooled'])} & "
    f"{'+' if x['delta_mean']>=0 else ''}{100*x['delta_mean']:.2f} \\\\" for k, x in v.items())

tex = r"""\section{Head-to-head against the published state of the art}
\label{sec:sota}

This section compares our methods with published rows on two public
benchmark datasets: on the
original AMP Scanner Vr.2 benchmark with its exact distributed splits, and
on the Feb2020 production benchmark with a 10-fold
cross-validation protocol on the same dataset.

\subsection{The Veltri 2018 benchmark}
Veltri, Kamath and Shehu's AMP Scanner Vr.2~\cite{ampscanner} is the most
cited deep-learning AMP recognizer, and Fu et al.'s ACEP~\cite{acep} reports
the strongest published numbers on the same data. The benchmark consists of
1{,}778 experimentally validated AMPs from the APD and 1{,}778 non-AMP
decoy peptides, identity-reduced and split once into 1{,}424 training, 708
tuning and 1{,}424 testing sequences; we downloaded the exact split files
from the authors' repository
(\texttt{github.com/dan-veltri/amp-scanner-v2}, \texttt{original-dataset/}).
ACEP's Table~2 is the canonical scoreboard for this test partition and we
reproduce it verbatim below, adding our own rows.

\begin{table}[h]\centering\small
\caption{Test-partition comparison on the Veltri 2018 benchmark. Published
rows are ACEP Table 2 verbatim; our stack is trained on the training
partition only, with stacking weights and the decision threshold tuned on
the tune partition (the role ACEP assigns it).}
\label{tab:sota}
\begin{tabular}{lccccc}
\hline
Method & SENS(\%) & SPEC(\%) & ACC(\%) & MCC & AUC(\%) \\
\hline
""" + tab1 + r"""
\hline
\end{tabular}
\end{table}

\subsection{Our protocol and per-model results}
We trained the full model zoo of this project on the benchmark's training
partition: three PepCNN seeds, two PepGNN seeds (shared-adjacency path) and
a random forest on a 426-dimensional descriptor vector (amino-acid and
dipeptide composition, net charge Eq.~\eqref{eq:netcharge}, mean and
standard deviation of Kyte--Doolittle hydropathy, Eisenberg moment
Eq.~\eqref{eq:hmoment}, Boman index, residue-class fractions). The stacking
logistic regression of Section~\ref{sec:stack-math} and the decision
threshold were tuned on the tune partition only.

\begin{table}[h]\centering\small
\caption{Base models and ensembles on the Veltri test partition.}
\label{tab:veltri-ours}
\begin{tabular}{lccc}
\hline
Model & ACC(\%) & MCC & AUC(\%) \\
\hline
""" + tab2 + r"""
\hline
\end{tabular}
\end{table}

The stacked ensemble reaches 91.19\% accuracy, MCC 0.825 and AUC 96.59\%:
close to the published AMP Scanner row (no paired significance test is available)
(91.29/0.8261/96.30) and behind ACEP (93.04/0.8610/97.78), which uses PSSM profiles. We
report a descriptive near-parity result, not a record.

\subsection{Testing against the Feb2020 production benchmark}
The production server model was retrained in February 2020 on an updated
dataset of 2{,}021 AMPs and 2{,}021 non-AMPs (downloaded from the authors'
site, \texttt{AMP\_Scan2\_Feb2020\_Dataset.zip}); its published 10-fold
cross-validation performance is SENS 90.6\%, SPEC 89.1\%, ACC 89.9\%, MCC
0.799, auROC 96.2\%. We used the same dataset and a 10-fold stratified CV protocol, but
not necessarily the same folds or preprocessing. Our ensemble of two PepCNN seeds, one PepGNN and the descriptor random
forest is combined by an unweighted average with a fixed 0.5 threshold, so no
tuning of any kind touches the held-out fold.

\begin{table}[h]\centering\small
\caption{10-fold CV on the complete Feb2020 dataset (4{,}042 sequences), vs
the published Feb2020 row.}
\label{tab:feb2020}
\begin{tabular}{cccccc}
\hline
Fold & SENS(\%) & SPEC(\%) & ACC(\%) & MCC & AUC(\%) \\
\hline
""" + fold_rows + r"""
\hline
\end{tabular}
\end{table}

\begin{table}[h]\centering\small
\caption{Verdict: published Feb2020 CV vs this work.}
\label{tab:feb2020-verdict}
\begin{tabular}{lcccc}
\hline
Metric & Published & Ours (mean) & Ours (pooled) & $\Delta$ (pp) \\
\hline
""" + verdict_tex + r"""
\hline
\end{tabular}
\end{table}

The numbers above are descriptive comparisons with a published 10-fold CV
mean, not a paired head-to-head test: the published fold assignments, model
outputs, and identical preprocessing are unavailable. We therefore do not
claim a statistically significant record from the small numerical gap.
"""
(ROOT/"paper/sota_benchmark.tex").write_text(tex)
print("sota_benchmark.tex written", len(tex))
