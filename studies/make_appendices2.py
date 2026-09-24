"""Generate additional paper appendices from committed result JSONs (study 18 support)."""
import json, pathlib
R = pathlib.Path('results'); P = pathlib.Path('paper')

def esc(s): return str(s).replace('_', r'\_')

# --- app_candidates_full.tex: all screen candidates with full audit columns ---
c = json.load(open(R/'study12_candidates.json'))
rows = []
for i, cand in enumerate(c['candidates'], 1):
    ph = cand['physicochemical']; nv = cand['novelty']
    rows.append(
        f"{i} & \\texttt{{{esc(cand.get('accession', cand['name']))}}} & {ph['length']} & "
        f"{cand['ensemble_score']:.4f} & {cand['cnn_score']:.4f} & {cand['gnn_score']:.4f} & "
        f"{nv['novelty_3mer_vs_APD6']:.3f} & {ph['net_charge_ph7.4']:.2f} & "
        f"{ph['kd_mean']:.2f} & {ph['hydrophobic_moment']:.2f} \\\\")
tab = "\n".join(rows)
seqs = []
for cand in c['candidates']:
    s = cand['sequence']
    wrapped = " ".join(s[j:j+10] for j in range(0, len(s), 10))
    org = cand.get('organism', 'annealed derivative of ' + cand.get('parent', '?'))
    seqs.append(f"\\item[\\texttt{{{esc(cand.get('accession', cand['name']))}}}] ({org}, {len(s)} aa) "
                f"\\texttt{{{wrapped}}}\n  {cand['falsifiable_prediction']}.")
seqblock = "\n".join(seqs)
(P/'app_candidates_full.tex').write_text(r"""
\section{Full discovery-screen candidate audit}\label{app:cand-full}
Table~\ref{tab:candfull} lists every candidate reported by the discovery
screen with its frozen-ensemble score, the CNN and GNN sub-scores, 3-mer
novelty against APD6 (Eq.~\eqref{eq:novelty}), and the physicochemical
audit values of Section~\ref{sec:physchem-math}. Scores are model
outputs under a confounded training set (Section~\ref{sec:labelconfound});
they are hypotheses for experimental testing, not measured activities.

\begin{table}[h]\centering\scriptsize
\caption{All discovery-screen candidates. Novelty is $1 - \max$ 3-mer
Jaccard vs APD6; $Q$ is net charge at pH 7.4; $\langle H\rangle$ is mean
Kyte--Doolittle hydropathy; $\mu_H$ is the Eisenberg moment.}
\label{tab:candfull}
\begin{tabular}{llcccccccc}
\hline
\# & Accession & Len & Ensemble & CNN & GNN & Novelty & $Q$ & $\langle H\rangle$ & $\mu_H$ \\
\hline
""" + tab + r"""
\hline
\end{tabular}
\end{table}

\subsection*{Sequences and falsifiable predictions}
Each candidate carries an explicit, testable prediction; a
broth-microdilution MIC assay against the named strains falsifies it.
\begin{description}
""" + seqblock + r"""
\end{description}
""")

# --- app_homology.tex: per-fold exact-decamer overlap ---
h = json.load(open(R/'study16_homology_audit.json'))
rows = "\n".join(f"{f['fold']} & {f['test_n']} & {f['shared_exact_10mer_n']} & {100*f['fraction']:.2f} \\\\" for f in h['folds'])
(P/'app_homology.tex').write_text(r"""
\section{Per-fold homology leakage audit}\label{app:homology}
Random stratified cross-validation does not control homology between the
training and held-out partitions. For each fold of the study-14 protocol we
count held-out sequences sharing at least one exact contiguous decamer
($k{=}10$) with any training sequence (the union-find criterion of
Definition~\ref{def:cluster} restricted to train--test pairs).

\begin{table}[h]\centering\small
\caption{Held-out sequences sharing an exact decamer with the training
partition, per fold of the Feb2020 10-fold CV.}
\label{tab:homologyfolds}
\begin{tabular}{cccc}
\hline
Fold & Held-out $n$ & Sharing a decamer & Fraction (\%) \\
\hline
""" + rows + r"""
\hline
\multicolumn{4}{l}{overall: """ + f"{100*h['overall_fraction']:.2f}" + r"""\% of all held-out sequences} \\
\end{tabular}
\end{table}

Roughly one in six held-out sequences is decamer-redundant with training,
so the random-CV numbers of Table~\ref{tab:feb2020} overstate performance
on genuinely remote sequences. A cluster-held-out variant of the same
protocol is the correct control and is recorded as future work; the
cluster-aware split of our own large dataset
(Section~\ref{sec:labelconfound}) shows the direction and size of the drop
such control induces.
""")

# --- app_strata_full.tex ---
s = json.load(open(R/'study17_review_strata.json'))
def strata_rows(d):
    out = []
    for model in ('PepCNN', 'PepGNN'):
        m = d[model]
        out.append(f"{model} & {m['auroc']:.4f} & {m['auprc']:.4f} & {m['mean_positive_score']:.4f} & {m['mean_negative_score']:.4f} \\\\")
    return "\n".join(out)
blocks = []
labels = {'all': 'All records (as trained)',
          'reviewed_only': 'Reviewed positives vs reviewed negatives',
          'unreviewed_positive_vs_reviewed_negative': 'Unreviewed positives vs reviewed negatives'}
for key in ('all', 'reviewed_only', 'unreviewed_positive_vs_reviewed_negative'):
    d = s[key]
    blocks.append(r"""
\begin{table}[h]\centering\small
\caption{Stratum: """ + labels[key] + f" ($n={d['n']:,}$, {d['n_pos']:,} positive, {d['n_neg']:,} negative)." + r"""}
\label{tab:strata-""" + key.replace('_','-') + r"""}
\begin{tabular}{lcccc}
\hline
Model & AUROC & AUPRC & Mean pos.\ score & Mean neg.\ score \\
\hline
""" + strata_rows(d) + r"""
\hline
\end{tabular}
\end{table}
""")
(P/'app_strata_full.tex').write_text(r"""
\section{Review-status stratification: complete tables}\label{app:strata}
The label-provenance audit of Section~\ref{sec:labelconfound} separates the
large-dataset test scores by UniProt review status. The tables below give
the complete numbers behind that analysis. The reviewed-positive stratum
is small (239 sequences) and its interval is correspondingly wide; the
unreviewed stratum dominates the training distribution and therefore the
headline score.
""" + "\n".join(blocks))

# --- app_sources.tex: full checksum manifest ---
b = json.load(open(R/'benchmark_sources.json'))
rows = []
for fn, meta in b['files'].items():
    rows.append(f"\\texttt{{{esc(fn)}}} & {meta['records']:,} & \\texttt{{{meta['sha256'][:16]}}} \\\\")
(P/'app_sources.tex').write_text(r"""
\section{Published-benchmark source manifest}\label{app:sources}
All published-benchmark files were downloaded from the original authors'
distribution points and pinned by SHA-256 before any experiment touched
them. Veltri 2018 splits: \texttt{github.com/dan-veltri/amp-scanner-v2},
\texttt{original-dataset/}. Feb2020 dataset:
\texttt{dveltri.com/ascan/v2}, \texttt{AMP\_Scan2\_Feb2020\_Dataset.zip}.
Full 64-hex hashes are in \texttt{results/benchmark\_sources.json}.

\begin{table}[h]\centering\small
\caption{Downloaded benchmark files: record counts and truncated SHA-256.}
\label{tab:manifest}
\begin{tabular}{lrl}
\hline
File & Records & SHA-256 (first 16) \\
\hline
""" + "\n".join(rows) + r"""
\hline
\end{tabular}
\end{table}
""")

# --- app_arch.tex: parameter tables (values verified against code) ---
(P/'app_arch.tex').write_text(r"""
\section{Architecture and parameter counts}\label{app:arch}
All parameter counts below are exact, verified against the implementation
(\texttt{pepdesign/models/}) with \texttt{count\_parameters}.

\subsection{PepCNN}
Each branch $b$ with kernel width $k_b$ applies two convolutional layers
with 32 channels and same-padding, followed by global max pooling over the
length axis. A convolution from $c_{\mathrm{in}}$ to $c_{\mathrm{out}}$
channels with kernel width $k$ costs
\begin{equation}
N_{\mathrm{conv}} = k\, c_{\mathrm{in}}\, c_{\mathrm{out}} + c_{\mathrm{out}}
\label{eq:convparams}
\end{equation}
parameters (weights plus biases). With $c_{\mathrm{in}} = 24$ input
channels (20 one-hot residues plus 4 physicochemical), widths
$k \in \{3,5,7\}$ and $c = 32$ hidden channels, each branch costs
$(24k \cdot 32 + 32) + (32k \cdot 32 + 32)$, and the head costs
$(96 \cdot 64 + 64) + (64 \cdot 1 + 1)$.

\begin{table}[h]\centering\small
\caption{PepCNN parameter budget (total 33{,}345).}
\label{tab:cnnparams}
\begin{tabular}{lccc}
\hline
Component & Shape & Parameters & Eq.~\eqref{eq:convparams} \\
\hline
Branch $k{=}3$, conv 1 & $32 \times 24 \times 3$ & 2{,}304 + 32 & $3\cdot24\cdot32+32$ \\
Branch $k{=}3$, conv 2 & $32 \times 32 \times 3$ & 3{,}072 + 32 & $3\cdot32\cdot32+32$ \\
Branch $k{=}5$, conv 1 & $32 \times 24 \times 5$ & 3{,}840 + 32 & $5\cdot24\cdot32+32$ \\
Branch $k{=}5$, conv 2 & $32 \times 32 \times 5$ & 5{,}120 + 32 & $5\cdot32\cdot32+32$ \\
Branch $k{=}7$, conv 1 & $32 \times 24 \times 7$ & 5{,}376 + 32 & $7\cdot24\cdot32+32$ \\
Branch $k{=}7$, conv 2 & $32 \times 32 \times 7$ & 7{,}168 + 32 & $7\cdot32\cdot32+32$ \\
Head linear 1 & $64 \times 96$ & 6{,}144 + 64 & \\
Head linear 2 & $1 \times 64$ & 64 + 1 & \\
\hline
\end{tabular}
\end{table}

The widest kernel sets the receptive field: after two $k{=}7$ layers each
residue logit sees $2(7-1)+1 = 13$ residues, a 13-mer context, close to
the 10-mer motif scale used by the homology audit.

\subsection{PepGNN}
A graph convolution $\mathrm{ReLU}(A_{\mathrm{norm}} H W + b)$ from $d$ to
$d'$ features costs $d d' + d'$ parameters. With 4 node features, three
layers of width 32 and the same 64-unit head, the model totals 4{,}449
parameters.

\begin{table}[h]\centering\small
\caption{PepGNN parameter budget (total 4{,}449).}
\label{tab:gnnparams}
\begin{tabular}{lccc}
\hline
Component & Shape & Parameters & Rule \\
\hline
GraphConv 1 & $4 \times 32$ & 128 + 32 & $dd' + d'$ \\
GraphConv 2 & $32 \times 32$ & 1{,}024 + 32 & $dd' + d'$ \\
GraphConv 3 & $32 \times 32$ & 1{,}024 + 32 & $dd' + d'$ \\
Head linear 1 & $64 \times 32$ & 2{,}048 + 64 & \\
Head linear 2 & $1 \times 64$ & 64 + 1 & \\
\hline
\end{tabular}
\end{table}

PepGNN is $7.5\times$ smaller than PepCNN, which partly explains its weaker
standalone performance (Tables \ref{tab:veltri-ours} and
\ref{tab:strata-all}): the chain-window adjacency carries strictly less
motif information than a 13-mer receptive field, and mean pooling dilutes
short local motifs that global max pooling preserves. Its marginal value
is diversity in the ensemble, not standalone accuracy.

\subsection{Descriptor random forest}
The classical baseline uses 400 trees (scikit-learn defaults otherwise)
over a 426-dimensional descriptor: 20 amino-acid fractions, 400 dipeptide
fractions, net charge \eqref{eq:netcharge}, mean and standard deviation of
Kyte--Doolittle hydropathy, the Eisenberg moment \eqref{eq:hmoment}, the
Boman index, and 3 residue-class fractions. The dipeptide block dominates
the dimension count, so the forest behaves like a regularized 2-mer model;
its competitive AUROC on the Veltri partition (96.19\%) is evidence that
much of the benchmark's signal is compositional rather than structural.
""")
print("appendices written")
