"""Generate paper tables (LaTeX) and figures from results/*.json."""
import json, pathlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parent.parent
RES, PAPER = ROOT / "results", ROOT / "paper"

s01 = json.load(open(RES / "study01_amp.json"))
s07 = json.load(open(RES / "study07_scaleup.json"))
s09 = json.load(open(RES / "study09_apd_validation.json"))
s08 = json.load(open(RES / "study08_discovery.json"))

def fmt(v): return f"{v:.3f}"

lines = []
# Table: pilot
lines.append("\\begin{table}[h]\\centering\\caption{Pilot benchmark (study01): 5-fold CV on 1{,}232 AMP vs 1{,}500 length-matched non-AMP, UniProt KW-0929, 8--60 aa. Brackets: bootstrap 95\\% CI; $p$: label-permutation.}\\label{tab:pilot}")
lines.append("\\begin{tabular}{lcccc}\\hline")
lines.append("Model & AUROC & 95\\% CI & AUPRC & $p$ \\\\ \\hline")
for m in ["PepCNN", "PepGNN", "logreg_dipeptide", "rf_dipeptide"]:
    v = s01[m]
    lines.append(f"{m.replace('_',' ')} & {fmt(v['auroc'])} & [{fmt(v['auroc_ci95'][0])}, {fmt(v['auroc_ci95'][1])}] & {fmt(v['auprc'])} & {v['permutation_p']:.4f} \\\\")
lines.append("\\hline\\end{tabular}\\end{table}")

# Table: scale-up
lines.append("\\begin{table}[h]\\centering\\caption{Scale-up benchmark (study07): held-out cluster-aware test split (4{,}386 sequences) from 24{,}628 AMP / 24{,}628 non-AMP; models trained on a 30{,}000-sequence subsample of the 40{,}776-sequence train split.}\\label{tab:scaleup}")
lines.append("\\begin{tabular}{lcccc}\\hline")
lines.append("Model & AUROC & 95\\% CI & AUPRC & $p$ \\\\ \\hline")
for m in ["PepCNN","PepGNN","logreg_dipeptide","rf_dipeptide"]:
    v = s07[m]
    lines.append(f"{m.replace('_',' ')} & {fmt(v['auroc'])} & [{fmt(v['auroc_ci95'][0])}, {fmt(v['auroc_ci95'][1])}] & {fmt(v['auprc'])} & {v['permutation_p']:.4f} \\\\")
lines.append("\\hline\\end{tabular}\\end{table}")

# Table: APD6 external
lines.append("\\begin{table}[h]\\centering\\caption{External validation (study09): APD6 2024 natural AMPs, decontaminated against train by shared 10-mers (%d of %d retained), scored against %d held-out test non-AMPs. Sensitivity at 5\\%% FPR: %.3f.}\\label{tab:apd}" % (s09["apd6_positives_decontaminated"], s09["apd6_positives_raw"], s09["negatives_test_split"], s09["apd6_sensitivity_at_test_fpr5"]))
lines.append("\\begin{tabular}{lccc}\\hline")
lines.append("Model & AUROC & 95\\% CI & AUPRC \\\\ \\hline")
for m, v in s09["models"].items():
    lines.append(f"{m.replace('_',' ')} & {fmt(v['auroc'])} & [{fmt(v['auroc_ci95'][0])}, {fmt(v['auroc_ci95'][1])}] & {fmt(v['auprc'])} \\\\")
lines.append("\\hline\\end{tabular}\\end{table}")

# Table: top designs (study07)
lines.append("\\begin{table}[h]\\centering\\caption{Top annealed designs (study07 ensemble; novelty = $1 - $ max 3-mer Jaccard vs AMP corpus). Scores are model confidences, not measured activities.}\\label{tab:designs}")
lines.append("\\begin{tabular}{lcc}\\hline Sequence & Score & Novelty \\\\ \\hline")
for d in s07["designs"][:10]:
    lines.append(f"\\texttt{{{d['sequence']}}} & {fmt(d['score'])} & {fmt(d['novelty_vs_corpus'])} \\\\")
lines.append("\\hline\\end{tabular}\\end{table}")

# Table: discovery top-15
lines.append("\\begin{table}[h]\\centering\\caption{Discovery screen (study08): top 15 of 852 reviewed uncharacterized 10--60 aa proteins by ensemble score. In-silico hypotheses only.}\\label{tab:screen}")
lines.append("\\begin{tabular}{llccc}\\hline Accession & Len & Ensemble & CNN & Novelty \\\\ \\hline")
for h in s08["top"][:15]:
    lines.append(f"{h['accession']} & {h['length']} & {fmt(h['ensemble_score'])} & {fmt(h['cnn_score'])} & {fmt(h['novelty_vs_amp_corpus'])} \\\\")
lines.append("\\hline\\end{tabular}\\end{table}")

(PAPER / "results_tables.tex").write_text("\n".join(lines) + "\n")
print("tables written")

# Figure: AUROC comparison bar chart
fig, ax = plt.subplots(figsize=(6.2, 3.4))
models = ["PepCNN", "PepGNN", "logreg", "rf"]
pilot = [s01["PepCNN"]["auroc"], s01["PepGNN"]["auroc"], s01["logreg_dipeptide"]["auroc"], s01["rf_dipeptide"]["auroc"]]
scale = [s07["PepCNN"]["auroc"], s07["PepGNN"]["auroc"], s07["logreg_dipeptide"]["auroc"], s07["rf_dipeptide"]["auroc"]]
apd = [s09["models"]["PepCNN"]["auroc"], s09["models"]["PepGNN"]["auroc"], float("nan"), float("nan")]
x = range(4); w = 0.27
ax.bar([i - w for i in x], pilot, w, label="Pilot CV (n=2,732)")
ax.bar(list(x), scale, w, label="Scale-up test (n=4,386)")
ax.bar([i + w for i in x], apd, w, label="APD6 external")
ax.set_xticks(list(x)); ax.set_xticklabels(["PepCNN", "PepGNN", "logreg", "RF"])
ax.set_ylim(0.7, 1.0); ax.set_ylabel("AUROC"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(PAPER / "fig_auroc.png", dpi=200)
print("figure written")
