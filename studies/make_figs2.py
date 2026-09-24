"""Figures for the benchmark/candidate sections."""
import json, pathlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parent.parent

# fig 1: SOTA comparison bar chart (ACC and MCC and AUC)
rows = [("iAMP-2L", 84.90, .6983, 84.90), ("CAMPr3-RF", 87.57, .7553, 93.63),
        ("gkmSVM", 89.46, .7895, 94.98), ("AMPScanner", 91.29, .8261, 96.30),
        ("ACEP", 93.04, .8610, 97.78)]
st = json.load(open(ROOT/"results/study11_veltri_stack.json"))["ensemble_stack_tuned"]
rows.append(("PepDesign stack", 100*st["acc"], st["mcc"], 100*st["auc"]))
names = [r[0] for r in rows]
x = range(len(names)); w = 0.28
fig, ax = plt.subplots(figsize=(7, 3.4))
ax.bar([i-w for i in x], [r[1] for r in rows], w, label="ACC %")
ax.bar([i for i in x], [100*r[2] for r in rows], w, label="MCC x100")
ax.bar([i+w for i in x], [r[3] for r in rows], w, label="AUC %")
ax.set_xticks(list(x)); ax.set_xticklabels(names, rotation=18, ha="right", fontsize=8)
ax.set_ylim(60, 100); ax.set_ylabel("%"); ax.legend(fontsize=8)
ax.set_title("Veltri 2018 test partition: published methods vs this work")
fig.tight_layout(); fig.savefig(ROOT/"paper/fig_sota.png", dpi=160); plt.close(fig)

# fig 2: study13 per-fold
s13 = json.load(open(ROOT/"results/study13_feb2020_cv.json"))
pub = s13["published_feb2020_cv"]
fs = s13["folds"]
x = [r["fold"] for r in fs]
fig, ax = plt.subplots(figsize=(7, 3.4))
ax.plot(x, [100*r["acc"] for r in fs], "o-", label="ACC % (ours)")
ax.plot(x, [100*r["auc"] for r in fs], "s-", label="AUC % (ours)")
ax.axhline(100*pub["acc"], color="tab:blue", ls="--", label="Feb2020 published ACC")
ax.axhline(100*pub["auc"], color="tab:orange", ls="--", label="Feb2020 published AUC")
ax.set_xlabel("fold"); ax.set_ylabel("%"); ax.set_xticks(x)
ax.set_ylim(80, 100); ax.legend(fontsize=8)
ax.set_title("10-fold CV on the AMP Scanner Feb2020 dataset")
fig.tight_layout(); fig.savefig(ROOT/"paper/fig_feb2020.png", dpi=160); plt.close(fig)
print("figs done")
