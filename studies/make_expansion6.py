"""Sixth expansion: segment table, environment, output schema, lessons."""
import json, pathlib
R = pathlib.Path('results'); P = pathlib.Path('paper')

# enrich case study with segment table
ab = json.load(open(R/'study15_signal_ablation.json'))
rows = []
for key, label in [('P0ACW4_precursor','Full precursor (1--57)'),
                   ('P0ACW4_signal_1_20','Signal segment (1--20)'),
                   ('P0ACW4_mature_after20','Post-signal segment (21--57)')]:
    v = ab[key]
    rows.append(f"{label} & {v['length']} & {v['cnn']:.4f} & {v['gnn']:.4f} & {v['ensemble']:.4f} \\\\")
seg_table = r"""
\begin{table}[h]\centering\small
\caption{Signal-peptide ablation for P0ACW4 with the frozen ensemble
(study 15; annotation: EcoliWiki ydcA gene record, signal residues 1--20).}
\label{tab:signalabl}
\begin{tabular}{lcccc}
\hline
Segment & Length & CNN & GNN & Ensemble \\
\hline
""" + "\n".join(rows) + r"""
\hline
\end{tabular}
\end{table}
"""
cs = P/'app_case.tex'
t = cs.read_text()
t = t.replace("\\paragraph{Verdict.}", seg_table + "\n\\paragraph{Verdict.}")
cs.write_text(t)

# app_env.tex
(P/'app_env.tex').write_text(r"""
\section{Compute environment}\label{app:env}
All experiments ran CPU-only in a 2-core Linux sandbox. Key library
versions: Python 3.11, PyTorch 2.14.0+cpu, scikit-learn 1.7.2, NumPy
2.2.6, SciPy 1.15.3. No GPU, no network at study time (fetches are
separate provenance-logged scripts). The full 10-fold study-14 run takes
about 15 minutes wall-clock (85--89 s per fold); the study-11 stack
training about 2 hours; the scale-up study-07 training is the multi-hour
step. Every reported number regenerates from the committed scripts and
pinned data on commodity hardware; nothing requires a cluster.

The hermetic test suite (27 tests) covers encoding invariants, model
shapes, the clustering union-find, the benchmark loaders against bundled
fixtures, the CLI, and a no-stubs audit that rejects placeholder
implementations. Tests run offline against cached fixtures and complete
in under two minutes.
""")

# app_schema.tex
(P/'app_schema.tex').write_text(r"""
\section{Screen output schema}\label{app:schema}
\texttt{pepdesign screen} writes one JSON document per run:

\begin{verbatim}
{
  "screened": <int>,               // sequences read from --input
  "top": [
    {
      "name": <str>,               // FASTA header, verbatim
      "sequence": <str>,           // cleaned, standard alphabet
      "length": <int>,
      "ensemble_score": <float>,   // mean of CNN and GNN probabilities
      "cnn_score": <float>,
      "gnn_score": <float>,
      "novelty_3mer": <float>,     // 1 - max 3-mer Jaccard vs reference
      "max_3mer_jaccard_vs_reference": <float>,
      "passes_novelty_floor": <bool>
    }, ...                         // sorted by ensemble_score desc
  ]
}
\end{verbatim}

Two contract points matter for downstream use. First, scores are frozen
ensemble outputs under the training distribution documented in
Section~\ref{sec:labelconfound}; consumers must treat them as rank orders,
not calibrated probabilities (Appendix~\ref{app:math3}). Second, novelty
is computed against the reference bundled with the model (APD6 by
default); screening against a different corpus requires re-running the
novelty audit, which \texttt{study12\_candidates.py} demonstrates against
both APD6 and the KW-0929 corpus.
""")

# lessons.tex (discussion addendum)
(P/'lessons.tex').write_text(r"""
\subsection{Lessons for AMP benchmark design}
Three features of this project generalize beyond the specific numbers.
First, \emph{review status is a label leak}: UniProt keyword annotation
without a \texttt{reviewed:true} filter produced a dataset whose positive
class is 88\% unreviewed while the negative class is 100\% reviewed, so a
model can profit from curation metadata rather than biology. Any UniProt
derived benchmark should report the review-status cross-tab alongside
accuracy. Second, \emph{length is a silent confound}: the reviewed
positives average 66 aa against 93--97 aa for the other groups
(Table~\ref{tab:compaudit}), and length-matching by construction, not just
by mean, is the cheap control that removes it. Third, \emph{random CV
flatters homology}: 15.98\% of held-out Feb2020 sequences share an exact
decamer with training (Appendix~\ref{app:homology}), so cross-validation
on this dataset estimates interpolation, not extrapolation to new
families. A published fold assignment or a cluster-held-out protocol is
the minimum requirement for a claim of the form ``method X beats method Y
on remote sequences.''
""")
print("expansion6 done")
