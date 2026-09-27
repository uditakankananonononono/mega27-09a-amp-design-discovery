"""Study 31 - Q10 encoding ablation (addendum 8, locked 2026-09-27 BEFORE
scoring). feb2020 2021+2021, identical study19 protocol: StratifiedKFold(10,
shuffle, random_state=42), CNN recipe verbatim (epochs 20, batch 64,
pos_weight 1.0, seeds 7/27). Arms: (a) full combined 24ch (comparator of
record), (b) one-hot-only 20ch, (c) physchem-only 4ch; CNN branches only per
lock (GNN/RF/HGB unchanged, not rerun). OOF MCC (thr 0.5) primary, AUC
secondary, paired folds, reported as-is; a drop is the expected outcome.

  PYTHONPATH=. python3 studies/study31_q10_encoding_ablation.py --arm full --folds 1-10
  PYTHONPATH=. python3 studies/study31_q10_encoding_ablation.py --arm onehot --folds 1-10
  PYTHONPATH=. python3 studies/study31_q10_encoding_ablation.py --arm physchem --folds 1-10
  PYTHONPATH=. python3 studies/study31_q10_encoding_ablation.py --finalize
"""
from __future__ import annotations
import argparse, json, math, pathlib, statistics as st
import numpy as np
import torch
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import matthews_corrcoef, roc_auc_score

from studies.study14_feb2020_full import load_fasta_seqs, train_binary, MAX_LEN
from pepdesign.models.cnn import PepCNN
from pepdesign.models.ensemble import predict_scores
from pepdesign.encoding import one_hot, physchem, combined_features

ROOT = pathlib.Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw/feb2020"
PART = ROOT / "results/study31_q10_partial"
FINAL = ROOT / "results/study31_q10_encoding_ablation.json"
FOLDS = 10
ARMS = {"full": (combined_features, 24), "onehot": (one_hot, 20), "physchem": (physchem, 4)}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", choices=list(ARMS))
    ap.add_argument("--folds", type=str, default="1-10")
    ap.add_argument("--finalize", action="store_true")
    args = ap.parse_args()
    amp = load_fasta_seqs(RAW/"AMPS_02182020.fasta")
    dec = load_fasta_seqs(RAW/"DECOYS_02182020.fasta")
    seqs = amp + dec
    y = np.array([1]*len(amp) + [0]*len(dec))
    assert (len(amp), len(dec)) == (2021, 2021)
    if args.finalize:
        out = {"study": "study31_q10_encoding_ablation",
               "protocol": "feb2020 2021+2021, study19 StratifiedKFold(10, shuffle, rs=42), CNN recipe verbatim, CNN branches only",
               "arms": {}}
        for arm in ARMS:
            rows = [json.load(open(p)) for p in sorted(PART.glob(f"{arm}_fold_*.json"),
                    key=lambda p: int(p.stem.split('_')[-1]))]
            assert len(rows) == FOLDS, (arm, len(rows))
            ent = {}
            for k in ("cnn7", "cnn27", "cnn_avg"):
                mccs = [r[k]["mcc"] for r in rows]; aucs = [r[k]["auc"] for r in rows]
                ent[k] = {"mcc_mean": round(st.mean(mccs), 4), "mcc_sd": round(st.stdev(mccs), 4),
                          "auc_mean": round(st.mean(aucs), 4), "fold_mcc": [round(x, 4) for x in mccs]}
            out["arms"][arm] = ent
        for k in ("cnn7", "cnn27", "cnn_avg"):
            full = out["arms"]["full"][k]["mcc_mean"]
            for arm in ("onehot", "physchem"):
                out["arms"][arm][k]["delta_mcc_vs_full"] = round(out["arms"][arm][k]["mcc_mean"] - full, 4)
        json.dump(out, open(FINAL, "w"), indent=1)
        print("FINAL:", json.dumps({arm: {k: v[k]["mcc_mean"] for k in ("cnn7","cnn27","cnn_avg")}
                                    for arm, v in out["arms"].items()}))
        return
    feat, dim = ARMS[args.arm]
    def inputs(ss, ml):
        return (torch.stack([torch.from_numpy(feat(x, ml)) for x in ss]),)
    skf = StratifiedKFold(n_splits=FOLDS, shuffle=True, random_state=42)
    splits = list(skf.split(seqs, y))
    PART.mkdir(exist_ok=True)
    lo, hi = (int(x) for x in args.folds.split("-"))
    for f in range(lo, hi+1):
        pf = PART / f"{args.arm}_fold_{f}.json"
        if pf.exists():
            print(f"{args.arm} fold {f}: exists, skip", flush=True); continue
        tri, tei = splits[f-1]
        tr_s = [seqs[i] for i in tri]; tr_y = y[tri]
        te_s = [seqs[i] for i in tei]; te_y = y[tei]
        row = {"arm": args.arm, "fold": f}
        ps = {}
        for name, seed in (("cnn7", 7), ("cnn27", 27)):
            cnn = train_binary(PepCNN(in_channels=dim), inputs, tr_s, tr_y, MAX_LEN,
                               epochs=20, batch_size=64, pos_weight=1.0, seed=seed)
            p = predict_scores(cnn, inputs, te_s, MAX_LEN)
            ps[name] = p
            pred = (p >= 0.5).astype(int)
            row[name] = {"mcc": float(matthews_corrcoef(te_y, pred)),
                         "auc": float(roc_auc_score(te_y, p))}
        pa = (ps["cnn7"] + ps["cnn27"]) / 2
        row["cnn_avg"] = {"mcc": float(matthews_corrcoef(te_y, (pa >= 0.5).astype(int))),
                          "auc": float(roc_auc_score(te_y, pa))}
        json.dump(row, open(pf, "w"))
        print(f"{args.arm} fold {f}: " + " ".join(f"{k} mcc {row[k]['mcc']:.4f}" for k in ("cnn7","cnn27","cnn_avg")), flush=True)

if __name__ == "__main__":
    main()
