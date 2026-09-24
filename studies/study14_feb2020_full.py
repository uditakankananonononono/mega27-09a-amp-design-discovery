"""Study 14 - full-record audit of the AMP Scanner Vr.2 Feb2020 production model
on its own benchmark (published 10-fold CV protocol).

Data: AMP_Scan2_Feb2020_Dataset.zip from dveltri.com/ascan/v2/news.html
(2,021 AMP + 2,021 non-AMP, the exact sequences the Feb2020 server model
was built on). Published Feb2020 10-fold CV numbers (news page):
  SENS 90.6% SPEC 89.1% ACC 89.9% MCC 0.799 auROC 96.2%
Protocol: same 10-fold stratified CV. Base learners: PepCNN x2 seeds,
PepGNN (shared adj) x1, RF on AAC+dipeptide+physicochemical descriptors.
Ensemble: unweighted average, threshold 0.5 (no tuning, fully nested).
Reported: mean +/- sd across folds for SENS/SPEC/ACC/MCC/AUROC, plus the
pooled out-of-fold metrics. Honest verdict vs the published row.
"""
from __future__ import annotations
import json, pathlib, time

import numpy as np
import torch
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import (roc_auc_score, accuracy_score, matthews_corrcoef)

from pepdesign.encoding import combined_features
from pepdesign.models.cnn import PepCNN
from pepdesign.models.gnn import PepGNN
from pepdesign.models.ensemble import (train_binary, predict_scores,
                                       train_binary_shared_adj,
                                       predict_scores_shared_adj)

ROOT = pathlib.Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw/feb2020"
MAX_LEN = 200
FOLDS = 10
CNN_SEEDS = (7, 27)

KD = {"A":1.8,"R":-4.5,"N":-3.5,"D":-3.5,"C":2.5,"Q":-3.5,"E":-3.5,"G":-0.4,
      "H":-3.2,"I":4.5,"L":3.8,"K":-3.9,"M":1.9,"F":2.8,"P":-1.6,"S":-0.8,
      "T":-0.7,"W":-0.9,"Y":-1.3,"V":4.2}
PKA = {"K":10.5,"R":12.4,"H":6.0,"D":3.9,"E":4.1,"C":8.3,"Y":10.1}
BOMAN = {"A":1.81,"R":-1.10,"N":0.76,"D":0.70,"C":1.09,"Q":0.67,"E":0.36,
         "G":1.14,"H":1.45,"I":2.59,"L":2.28,"K":-0.99,"M":1.84,"F":2.70,
         "P":-0.99,"S":0.54,"T":0.96,"W":3.07,"Y":2.76,"V":2.59}
AA = "ACDEFGHIKLMNPQRSTVWY"

def net_charge(seq, ph=7.4):
    pos = sum(10**PKA[a]/(10**ph+10**PKA[a]) for a in seq if a in ("K","R","H"))
    neg = sum(10**ph/(10**ph+10**PKA[a]) for a in seq if a in ("D","E","C","Y"))
    return pos - neg + 10**7.0/(10**ph+10**7.0) - 10**ph/(10**ph+10**2.34)

def hydrophobic_moment(seq, angle=100.0, window=11):
    if len(seq) < window: window = len(seq)
    rad = np.deg2rad(angle); best = 0.0
    for i in range(len(seq)-window+1):
        w = seq[i:i+window]
        x = sum(KD.get(a,0.0)*np.cos(j*rad) for j,a in enumerate(w))
        y = sum(KD.get(a,0.0)*np.sin(j*rad) for j,a in enumerate(w))
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
    kd = [KD.get(a,0.0) for a in seq]
    v += [n, net_charge(seq), float(np.mean(kd)), float(np.std(kd)),
          hydrophobic_moment(seq), float(np.mean([BOMAN.get(a,0.0) for a in seq]))]
    return v

def load_fasta_seqs(path):
    seqs, cur = [], ""
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if cur: seqs.append(cur)
                cur = ""
            elif line:
                cur += line
    if cur: seqs.append(cur)
    return seqs

def cnn_inputs(seqs, max_len):
    return (torch.stack([torch.from_numpy(combined_features(x, max_len)) for x in seqs]),)

def main():
    amp = load_fasta_seqs(RAW/"AMPS_02182020.fasta")
    dec = load_fasta_seqs(RAW/"DECOYS_02182020.fasta")
    seqs = amp + dec
    y = np.array([1]*len(amp) + [0]*len(dec))
    print(f"dataset: {len(amp)} AMP + {len(dec)} decoy = {len(seqs)}", flush=True)
    assert (len(amp),len(dec)) == (2021,2021)
    assert max(map(len,seqs)) <= MAX_LEN

    published = {"sens":0.906,"spec":0.891,"acc":0.899,"mcc":0.799,"auc":0.962}
    results = {"benchmark": "AMP Scanner Vr.2 Feb2020 dataset (dveltri.com)",
               "n_amp": len(amp), "n_decoy": len(dec),
               "protocol": "10-fold stratified CV, avg ensemble (2x PepCNN + PepGNN + RF-rich), thr 0.5, all 4042 original sequences retained (X treated as unknown)",
               "published_feb2020_cv": published}

    oof_sum = np.zeros(len(seqs)); oof_cnt = np.zeros(len(seqs))
    fold_rows = []
    skf = StratifiedKFold(n_splits=FOLDS, shuffle=True, random_state=42)
    for fold, (tri, tei) in enumerate(skf.split(seqs, y), 1):
        t0 = time.time()
        tr_s = [seqs[i] for i in tri]; tr_y = y[tri]
        te_s = [seqs[i] for i in tei]; te_y = y[tei]
        scores = []
        for seed in CNN_SEEDS:
            cnn = train_binary(PepCNN(), cnn_inputs, tr_s, tr_y, MAX_LEN,
                               epochs=20, batch_size=64, pos_weight=1.0, seed=seed)
            scores.append(predict_scores(cnn, cnn_inputs, te_s, MAX_LEN))
        gnn = train_binary_shared_adj(PepGNN(), tr_s, tr_y, MAX_LEN, epochs=10,
                                      batch_size=64, pos_weight=1.0, seed=11)
        scores.append(predict_scores_shared_adj(gnn, te_s, MAX_LEN))
        rf = RandomForestClassifier(n_estimators=300, n_jobs=2, random_state=42)
        rf.fit(np.array([rich_features(s) for s in tr_s]), tr_y)
        scores.append(rf.predict_proba(np.array([rich_features(s) for s in te_s]))[:,1])
        ens = np.mean(scores, axis=0)
        oof_sum[tei] += ens; oof_cnt[tei] += 1
        pred = (ens >= 0.5).astype(int)
        tp=int(((pred==1)&(te_y==1)).sum()); fn=int(((pred==0)&(te_y==1)).sum())
        tn=int(((pred==0)&(te_y==0)).sum()); fp=int(((pred==1)&(te_y==0)).sum())
        row = {"fold": fold, "sens": tp/(tp+fn), "spec": tn/(tn+fp),
               "acc": accuracy_score(te_y, pred),
               "mcc": float(matthews_corrcoef(te_y, pred)),
               "auc": roc_auc_score(te_y, ens),
               "seconds": round(time.time()-t0,1)}
        fold_rows.append(row)
        print(f"fold {fold}: acc {row['acc']:.3f} mcc {row['mcc']:.3f} auc {row['auc']:.3f} ({row['seconds']}s)", flush=True)

    oof = oof_sum / oof_cnt
    pred = (oof >= 0.5).astype(int)
    tp=int(((pred==1)&(y==1)).sum()); fn=int(((pred==0)&(y==1)).sum())
    tn=int(((pred==0)&(y==0)).sum()); fp=int(((pred==1)&(y==0)).sum())
    pooled = {"sens": tp/(tp+fn), "spec": tn/(tn+fp),
              "acc": accuracy_score(y, pred), "mcc": float(matthews_corrcoef(y, pred)),
              "auc": roc_auc_score(y, oof)}
    mean_sd = {m: [float(np.mean([r[m] for r in fold_rows])),
                   float(np.std([r[m] for r in fold_rows]))]
               for m in ("sens","spec","acc","mcc","auc")}
    results["folds"] = fold_rows
    results["mean_sd"] = mean_sd
    results["pooled_oof"] = pooled
    results["verdict"] = {k: {"published": published[k], "ours_mean": mean_sd[k][0],
                              "ours_pooled": pooled[k],
                              "beat_on_mean": mean_sd[k][0] > published[k],
                              "delta_mean": mean_sd[k][0] - published[k]}
                          for k in published}
    print(json.dumps(results["verdict"], indent=2), flush=True)
    out = ROOT/"results/study14_feb2020_full.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"wrote {out}", flush=True)

if __name__ == "__main__":
    main()
