"""Fourth expansion: encoding reference, designer diagnostics, case study."""
import json, pathlib
R = pathlib.Path('results'); P = pathlib.Path('paper')

# app_encoding.tex with exact per-residue values from code
import sys; sys.path.insert(0, '.')
from pepdesign.encoding import AA, HYDROPHOBICITY, CHARGE, VOLUME, AROMATIC
rows = []
for a in AA:
    rows.append(f"{a} & {HYDROPHOBICITY[a]:+.1f} & {HYDROPHOBICITY[a]/4.5:+.3f} & "
                f"{CHARGE[a]:+.1f} & {VOLUME[a]:.1f} & {VOLUME[a]/230.0:.3f} & "
                f"{'1' if a in AROMATIC else '0'} \\\\")
(P/'app_encoding.tex').write_text(r"""
\section{Encoding reference}\label{app:encoding}
The CNN consumes 24 channels per residue: the 20-indicator one-hot vector
followed by the four normalized physicochemical channels of
Table~\ref{tab:encoding}. The GNN consumes only the four physicochemical
channels as node features. Normalization keeps every channel in
$[-1, 1]$: hydropathy spans $[-4.5, 4.5]$, volume spans $[60.1, 227.8]$.

\begin{table}[h]\centering\small
\caption{Per-residue feature values, exactly as implemented in
\texttt{pepdesign/encoding.py}. $h$ = Kyte--Doolittle hydropathy,
$q$ = formal charge at physiological pH, $v$ = residue volume in \AA$^3$.}
\label{tab:encoding}
\begin{tabular}{lcccccc}
\hline
Residue & $h$ & $h/4.5$ & $q$ & $v$ & $v/230$ & aromatic \\
\hline
""" + "\n".join(rows) + r"""
\hline
\end{tabular}
\end{table}

Sequences are truncated or zero-padded to the study's maximum length $L$
(pilot 60, large studies 150, study 14: 200); non-standard residues are
dropped by \texttt{clean\_sequence}, and study 14 retains records
containing \texttt{X} by treating them as unknowns removed from the
encoding, a choice recorded because the published Feb2020 dataset contains
such records. Padding rows are all-zero in every channel, so the global
max pool of the CNN and the mean pool of the GNN see them as neutral
except where a learned negative bias activates; the head therefore never
receives an explicit length signal, and length enters only through the
pooled features.

The GNN adjacency is the chain window with self loops: $A_{ij} = 1$ iff
$|i-j| \le w$ or $i = j$ (window $w=2$ unless noted), row-normalized as
$D^{-1}A$ with $D_{ii} = \sum_j A_{ij}$. The shared-adjacency variant
broadcasts one $(L,L)$ matrix across the batch; the docstring in
\texttt{pepdesign/models/gnn.py} records the memory argument
($O(nL^2)$ per-sample adjacency is about 4.5 GB at $L{=}150$, $n{=}50$k,
which does not fit the sandbox).
""")

# app_designer.tex
s = json.load(open(R/'study07_scaleup.json'))
d = s['designs']
rows = []
for i, x in enumerate(d, 1):
    seq = x['sequence']
    rows.append(f"{i} & \\texttt{{{seq}}} & {x['score']:.4f} & {x['novelty_vs_corpus']:.3f} \\\\")
(P/'app_designer.tex').write_text(r"""
\section{Designer diagnostics}\label{app:designer}
The annealing designer of the methods section optimizes the ensemble score
with a novelty floor. Table~\ref{tab:alldesigns} lists the top-10 anneals;
Figure~\ref{fig:anneal} shows the score--novelty plane. The anneals
cluster at novelty 0.89--0.96 against the training corpus: the novelty
floor is active, and the optimizer is not simply memorizing a training
motif. All designs are short (25--57 aa), cationic and hydrophobic by
construction of the scoring function; none is a claimed AMP, and each
carries the same training-label confound as the screen candidates
(Section~\ref{sec:labelconfound}).

\begin{figure}[h]\centering
\includegraphics[width=0.65\textwidth]{fig_anneal.png}
\caption{Annealed designs in the score--novelty plane.}
\label{fig:anneal}
\end{figure}
""")

# app_case.tex: end-to-end P0ACW4 case study
c = json.load(open(R/'study12_candidates.json'))['candidates'][0]
ab = json.load(open(R/'study15_signal_ablation.json'))
(P/'app_case.tex').write_text(r"""
\section{Case study: one candidate end to end}\label{app:case}
This appendix traces the top screen hit through every stage, collecting
numbers reported piecemeal above into one falsifiable narrative.

\paragraph{Stage 1: pool.} The discovery pool is 852 reviewed,
uncharacterized UniProt proteins of 10--60 residues, fetched with
provenance logging. P0ACW4 (YdcA, \emph{E. coli} K-12, 57 aa) enters as an
uncharacterized protein; the screen makes no use of any AMP annotation.

\paragraph{Stage 2: scoring.} The frozen study-07 ensemble scores it
0.9082 (CNN 0.9992, GNN 0.8171), the top unique score in the pool. The
duplicate-accession pair P0ACW4/P0ACW5 carries an identical 57-residue
sequence in two species records and is reported as such.

\paragraph{Stage 3: novelty audit.} Maximum 3-mer Jaccard against APD6
2024 is 0.062 (novelty 0.938) and against the KW-0929 corpus 0.055; the
longest common substring with any APD6 entry is 6 residues. Under
Definition of Eq.~\eqref{eq:novelty} this is weak local similarity only.

\paragraph{Stage 4: physicochemical audit.} Length 57, net charge
$+4.87$ at pH 7.4 (Eq.~\eqref{eq:netcharge}), mean hydropathy $+0.08$,
Eisenberg moment $1.18$ (Eq.~\eqref{eq:hmoment}), Boman index $1.09$:
cationic and roughly hydrophobicity-neutral, a plausible but not extreme
AMP-like profile.

\paragraph{Stage 5: segment ablation.} YdcA carries an annotated 20-aa
N-terminal signal peptide (EcoCyc/EcoliWiki gene record). Scoring the
segments separately with the same frozen ensemble: full precursor 0.9082;
precursor minus signal segment 0.6299. The model's confidence is carried
substantially by the signal segment, which is expected: signal peptides
are hydrophobic, and the confounded training set rewards hydrophobic
cationic stretches (Section~\ref{sec:labelconfound}).

\paragraph{Verdict.} P0ACW4 is a named, quantified, falsifiable
\emph{hypothesis}: broth-microdilution MIC against \emph{E. coli} K-12 and
\emph{S. aureus} ATCC 25923 decides it. The ablation is precisely the kind
of negative control that keeps this report honest: the high score partly
reflects a trafficking signal, not necessarily antimicrobial chemistry.
""")
print("expansion4 done")
