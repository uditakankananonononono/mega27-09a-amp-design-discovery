"""Study 30 - Q9 negative controls through the discovery screen (addendum 8,
locked 2026-09-27 BEFORE scoring). Arms: (a) scrambles - each of the 10 named
InstiAMP-21 sequences shuffled 10x rng(11) (length+composition preserved);
(b) known non-AMP - 50 reviewed non-AMP proteins from the Q1 reviewed-negative
pool (data/raw/q1_matched_negatives.fasta), rng(11) sample. Screen: the locked
study21 rerun ensemble recipe (same 5 branches, same seeds/epochs, retrained
on full Feb2020 identically to study21_rescreen.py - training code copied
verbatim from that locked path). Locked expectation, not a threshold: controls
should score low. Full distributions reported; every p >= 0.5 case listed
individually, as-is; high scorers never dropped.
"""
from __future__ import annotations
import json, pathlib
import numpy as np

from studies.study14_feb2020_full import (rich_features,
    load_fasta_seqs, train_binary, train_binary_shared_adj, MAX_LEN)
from studies.study19_feb2020_v2 import rich_features_v2
from pepdesign.models.cnn import PepCNN
from pepdesign.models.gnn import PepGNN
from pepdesign.models.ensemble import predict_scores, predict_scores_shared_adj
from pepdesign.encoding import combined_features
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier

ROOT = pathlib.Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw/feb2020"
WEIGHTS = {"cnn7": .20, "cnn27": .20, "gnn": .15, "rf": .20, "hgb": .25}
OUT = ROOT / "results/study30_q9_negative_controls.json"

def cnn_inputs(seqs, max_len):
    return (torch_stack(seqs, max_len),)

import torch
def torch_stack(seqs, max_len):
    return torch.stack([torch.from_numpy(combined_features(x, max_len)) for x in seqs])

def main():
    amp = load_fasta_seqs(RAW/"AMPS_02182020.fasta")
    dec = load_fasta_seqs(RAW/"DECOYS_02182020.fasta")
    seqs = amp + dec
    y = np.array([1]*len(amp) + [0]*len(dec))
    assert (len(amp), len(dec)) == (2021, 2021)
    print("training study21-recipe ensemble on full Feb2020 ...", flush=True)
    branches = {}
    for name, seed in (("cnn7", 7), ("cnn27", 27)):
        branches[name] = train_binary(PepCNN(), cnn_inputs, seqs, y, MAX_LEN,
                                      epochs=20, batch_size=64, pos_weight=1.0, seed=seed)
        print(f"  {name} trained", flush=True)
    branches["gnn"] = train_binary_shared_adj(PepGNN(), seqs, y, MAX_LEN, epochs=10,
                                              batch_size=64, pos_weight=1.0, seed=11)
    print("  gnn trained", flush=True)
    X = np.array([rich_features(s) for s in seqs])
    rf = RandomForestClassifier(n_estimators=300, n_jobs=2, random_state=42).fit(X, y)
    X2 = np.array([rich_features_v2(s) for s in seqs])
    hgb = HistGradientBoostingClassifier(random_state=42).fit(X2, y)
    print("  rf+hgb trained", flush=True)

    def ens(s):
        c7 = predict_scores(branches["cnn7"], cnn_inputs, [s], MAX_LEN)[0]
        c27 = predict_scores(branches["cnn27"], cnn_inputs, [s], MAX_LEN)[0]
        g = predict_scores_shared_adj(branches["gnn"], [s], MAX_LEN)[0]
        r = rf.predict_proba(np.array([rich_features(s)]))[:, 1][0]
        h = hgb.predict_proba(np.array([rich_features_v2(s)]))[:, 1][0]
        return float(.20*c7 + .20*c27 + .15*g + .20*r + .25*h)

    rescreen = json.load(open(ROOT/"results/study21_rescreen.json"))
    named = [c["sequence"] for c in rescreen["named_candidates"]]
    rng = np.random.RandomState(11)
    scrambles = []
    for i, s in enumerate(named, 1):
        for k in range(10):
            scrambles.append({"name": f"scramble-21-{i}-{k+1}",
                              "sequence": "".join(rng.permutation(list(s)))})
    negsrc = load_fasta_seqs(ROOT/"data/raw/q1_matched_negatives.fasta")
    negsrc = [s for s in negsrc if 10 <= len(s) <= 150]
    pick = rng.choice(len(negsrc), size=50, replace=False)
    nonamp = [{"name": f"nonamp-{k+1}", "sequence": negsrc[j]} for k, j in enumerate(pick)]
    print(f"scoring {len(scrambles)} scrambles + {len(nonamp)} known non-AMP ...", flush=True)

    def score_arm(items):
        ps = [ens(it["sequence"]) for it in items]
        for it, p in zip(items, ps):
            it["ensemble_p"] = round(p, 4)
        arr = np.array(ps)
        return {"n": len(ps), "median": round(float(np.median(arr)), 4),
                "q1": round(float(np.percentile(arr, 25)), 4),
                "q3": round(float(np.percentile(arr, 75)), 4),
                "max": round(float(arr.max()), 4),
                "n_ge_0.5": int((arr >= 0.5).sum()),
                "ge_0.5_cases": [{"name": it["name"], "p": it["ensemble_p"]}
                                 for it in items if it["ensemble_p"] >= 0.5],
                "items": items}
    out = {"study": "study30_q9_negative_controls",
           "screen": "locked study21 rerun ensemble recipe (5 branches, same seeds/epochs, full Feb2020 retrain)",
           "arms": {"scrambles": score_arm(scrambles), "known_non_amp": score_arm(nonamp)}}
    json.dump(out, open(OUT, "w"), indent=1)
    for arm, r in out["arms"].items():
        print(arm, {k: r[k] for k in ("n","median","q1","q3","max","n_ge_0.5")}, flush=True)
        for c in r["ge_0.5_cases"]:
            print("  HIGH:", c, flush=True)

if __name__ == "__main__":
    main()
