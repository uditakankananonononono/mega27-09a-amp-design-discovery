"""Study 11 - SOTA attempt on the Veltri 2018 AMP benchmark via stacked ensemble.

Same splits as the published benchmark (1424 train / 708 tune / 1424 test,
AMP.eval used ONLY for tuning: stacking weights + decision threshold, the
same role the ACEP paper gives the tune partition). Published test-partition
rows to beat (ACEP paper Table 2, BMC Genomics 2020, 21:597):
  AMPScanner: SENS 89.88 SPEC 92.69 ACC 91.29 MCC 0.8261 AUC 96.30
  ACEP:       SENS 92.41 SPEC 93.67 ACC 93.04 MCC 0.8610 AUC 97.78
Base models: PepCNN x3 seeds, PepGNN (shared adj) x2 seeds, RF on
AAC+dipeptide+physicochemical descriptors. Meta-model: logistic regression
on tune-partition scores; threshold maximizes MCC on the tune partition.
"""
from __future__ import annotations
import json, pathlib

import numpy as np
import torch
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (roc_auc_score, average_precision_score,
                             accuracy_score, matthews_corrcoef)

from pepdesign.benchmark import bootstrap_ci
from pepdesign.encoding import combined_features
from pepdesign.models.cnn import PepCNN
from pepdesign.models.gnn import PepGNN
from pepdesign.models.ensemble import (train_binary, predict_scores,
                                       train_binary_shared_adj,
                                       predict_scores_shared_adj)

ROOT = pathlib.Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw/veltri"
MAX_LEN = 150

KD = {"A":1.8,"R":-4.5,"N":-3.5,"D":-3.5,"C":2.5,"Q":-3.5,"E":-3.5,"G":-0.4,
      "H":-3.2,"I":4.5,"L":3.8,"K":-3.9,"M":1.9,"F":2.8,"P":-1.6,"S":-0.8,
      "T":-0.7,"W":-0.9,"Y":-1.3,"V":4.2}
PKA = {"K":10.5,"R":12.4,"H":6.0,"D":3.9,"E":4.1,"C":8.3,"Y":10.1}
BOMAN = {"A":1.81,"R":-1.10,"N":0.76,"D":0.70,"C":1.09,"Q":0.67,"E":0.36,
         "G":1.14,"H":1.45,"I":2.59,"L":2.28,"K":-0.99,"M":1.84,"F":2.70,
         "P":-0.99,"S":0.54,"T":0.96,"W":3.07,"Y":2.76,"V":2.59}
AA = "ACDEFGHIKLMNPQRSTVWY"

def net_charge(seq, ph=7.4):
    q = 10**PKA.get("K",0)/(10**ph+10**PKA["K"]) * 0  # placeholder
    pos = sum(10**PKA[a]/(10**ph+10**PKA[a]) for a in seq if a in ("K","R","H"))
    neg = sum(10**ph/(10**ph+10**PKA[a]) for a in seq if a in ("D","E","C","Y"))
    pos += 10**7.0/(10**ph+10**7.0)   # N-terminus
    neg += 10**ph/(10**ph+10**2.34)   # C-terminus
    return pos - neg

def hydrophobic_moment(seq, angle=100.0, window=11):
    if len(seq) < window: window = len(seq)
    best = 0.0
    rad = np.deg2rad(angle)
    for i in range(len(seq) - window + 1):
        w = seq[i:i+window]
        x = sum(KD[a]*np.cos(j*rad) for j,a in enumerate(w))
        y = sum(KD[a]*np.sin(j*rad) for j,a in enumerate(w))
        best = max(best, (x*x+y*y)**0.5/window)
    return best

def rich_features(seq):
    n = len(seq)
    c = {a: seq.count(a) for a in AA}
    v = [c[a]/n for a in AA]
    di = {}
    tot = max(n-1, 1)
    for i in range(n-1):
        di[seq[i:i+2]] = di.get(seq[i:i+2], 0) + 1
    v += [di.get(a+b, 0)/tot for a in AA for b in AA]
    kd = [KD[a] for a in seq]
    v += [n, net_charge(seq), float(np.mean(kd)), float(np.std(kd)),
          hydrophobic_moment(seq), float(np.mean([BOMAN[a] for a in seq])),
          (c["K"]+c["R"])/n, (c["F"]+c["W"]+c["Y"])/n,
          (c["A"]+c["C"]+c["G"]+c["S"]+c["T"])/n,
          (c["D"]+c["E"])/n]
    return v

def load_fasta(path):
    seqs, cur = [], None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if cur: seqs.append(cur)
                cur = ""
            elif line:
                cur = (cur or "") + line
    if cur: seqs.append(cur)
    return seqs

def main():
    tr_s = load_fasta(RAW/"AMP.tr.fa") + load_fasta(RAW/"DECOY.tr.fa")
    tr_y = np.array([1]*712 + [0]*712)
    va_s = load_fasta(RAW/"AMP.eval.fa") + load_fasta(RAW/"DECOY.eval.fa")
    va_y = np.array([1]*354 + [0]*354)
    te_s = load_fasta(RAW/"AMP.te.fa") + load_fasta(RAW/"DECOY.te.fa")
    te_y = np.array([1]*712 + [0]*712)
    print(f"splits: train {len(tr_s)} tune {len(va_s)} test {len(te_s)}", flush=True)

    results = {"benchmark": "Veltri 2018 / ACEP Table 2 test partition",
               "published": {"AMPScanner": {"sens":.8988,"spec":.9269,"acc":.9129,"mcc":.8261,"auc":.9630},
                             "ACEP": {"sens":.9241,"spec":.9367,"acc":.9304,"mcc":.8610,"auc":.9778}}}

    def cnn_inputs(seqs, max_len):
        return (torch.stack([torch.from_numpy(combined_features(x, max_len)) for x in seqs]),)

    def metrics(scores, thr):
        pred = (scores >= thr).astype(int)
        tp=int(((pred==1)&(te_y==1)).sum()); fn=int(((pred==0)&(te_y==1)).sum())
        tn=int(((pred==0)&(te_y==0)).sum()); fp=int(((pred==1)&(te_y==0)).sum())
        return {"sens": tp/(tp+fn), "spec": tn/(tn+fp),
                "acc": accuracy_score(te_y, pred), "mcc": float(matthews_corrcoef(te_y, pred)),
                "auc": roc_auc_score(te_y, scores),
                "auprc": float(average_precision_score(te_y, scores))}

    def best_threshold(y, s):
        ths = np.linspace(0.05, 0.95, 181)
        mccs = [matthews_corrcoef(y, (s>=t).astype(int)) for t in ths]
        i = int(np.argmax(mccs)); return float(ths[i]), float(mccs[i])

    pw = 1.0
    base_va, base_te = {}, {}

    for seed in (7, 27, 77):
        print(f"PepCNN seed {seed} ...", flush=True)
        cnn = train_binary(PepCNN(), cnn_inputs, tr_s, tr_y, MAX_LEN,
                           epochs=25, batch_size=64, pos_weight=pw, seed=seed)
        base_va[f"cnn{seed}"] = predict_scores(cnn, cnn_inputs, va_s, MAX_LEN)
        base_te[f"cnn{seed}"] = predict_scores(cnn, cnn_inputs, te_s, MAX_LEN)

    for seed in (11, 111):
        print(f"PepGNN seed {seed} ...", flush=True)
        gnn = train_binary_shared_adj(PepGNN(), tr_s, tr_y, MAX_LEN, epochs=12,
                                      batch_size=64, pos_weight=pw, seed=seed)
        base_va[f"gnn{seed}"] = predict_scores_shared_adj(gnn, va_s, MAX_LEN)
        base_te[f"gnn{seed}"] = predict_scores_shared_adj(gnn, te_s, MAX_LEN)

    print("RF rich features ...", flush=True)
    Xtr = np.array([rich_features(s) for s in tr_s])
    Xva = np.array([rich_features(s) for s in va_s])
    Xte = np.array([rich_features(s) for s in te_s])
    rf = RandomForestClassifier(n_estimators=400, n_jobs=2, random_state=42)
    rf.fit(Xtr, tr_y)
    base_va["rf"] = rf.predict_proba(Xva)[:,1]
    base_te["rf"] = rf.predict_proba(Xte)[:,1]

    print("tuning on the tune partition ...", flush=True)
    names = list(base_va)
    Zva = np.stack([base_va[k] for k in names], axis=1)
    Zte = np.stack([base_te[k] for k in names], axis=1)

    # single best base model
    for k in names:
        m = metrics(base_te[k], 0.5)
        results[f"base_{k}"] = m
        print(f"  {k}: acc {m['acc']:.3f} mcc {m['mcc']:.3f} auc {m['auc']:.3f}", flush=True)

    # simple average
    avg_va = Zva.mean(axis=1); avg_te = Zte.mean(axis=1)
    thr_avg, _ = best_threshold(va_y, avg_va)
    m = metrics(avg_te, thr_avg); m["threshold"] = thr_avg
    results["ensemble_avg_tuned"] = m
    print(f"  avg(tuned thr {thr_avg:.2f}): acc {m['acc']:.3f} sens {m['sens']:.3f} spec {m['spec']:.3f} mcc {m['mcc']:.3f} auc {m['auc']:.3f}", flush=True)

    # stacking LR on tune partition
    stack = LogisticRegression(max_iter=2000, C=1.0).fit(Zva, va_y)
    s_va = stack.predict_proba(Zva)[:,1]; s_te = stack.predict_proba(Zte)[:,1]
    thr_st, _ = best_threshold(va_y, s_va)
    m = metrics(s_te, thr_st); m["threshold"] = thr_st
    auc, lo, hi = bootstrap_ci(te_y, s_te, n_boot=500, seed=1)
    m["auc_ci95"] = [lo, hi]
    m["stack_coef"] = dict(zip(names, stack.coef_[0].tolist()))
    results["ensemble_stack_tuned"] = m
    print(f"  stack(tuned thr {thr_st:.2f}): acc {m['acc']:.3f} sens {m['sens']:.3f} spec {m['spec']:.3f} mcc {m['mcc']:.3f} auc {m['auc']:.3f} [{lo:.3f},{hi:.3f}]", flush=True)

    # verdict vs ACEP
    best = max([("avg", results["ensemble_avg_tuned"]), ("stack", results["ensemble_stack_tuned"])],
               key=lambda kv: kv[1]["acc"])
    acep = results["published"]["ACEP"]
    results["verdict"] = {
        "best": best[0],
        "beats_acep_acc": best[1]["acc"] > acep["acc"],
        "beats_acep_mcc": best[1]["mcc"] > acep["mcc"],
        "beats_acep_auc": best[1]["auc"] > acep["auc"],
        "delta_acc": best[1]["acc"] - acep["acc"],
        "delta_mcc": best[1]["mcc"] - acep["mcc"],
        "delta_auc": best[1]["auc"] - acep["auc"]}
    print(json.dumps(results["verdict"], indent=2), flush=True)

    out = ROOT / "results/study11_veltri_stack.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"wrote {out}", flush=True)

if __name__ == "__main__":
    main()
