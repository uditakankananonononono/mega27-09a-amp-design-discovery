"""Study 21 - discovery re-screen with the study19 locked ensemble (addendum 1
discovery arm + addendum 2 A2.2 tiered novelty audit; both locked pre-run).
Trains the identical 5 branches on the FULL Feb2020 dataset (same configs and
seeds as the study19 folds; full-data training is forced by the locked arm
design - the ensemble that passed gates is a CV artifact, the screen needs a
final model), then:
 1. scores the study08 uncharacterized UniProt pool (852 seqs) and study12's
    9 named candidates;
 2. tiered novelty audit (locked A2.2): 3-mer Jaccard + LCS vs APD6 2024,
    DRAMP natural (pinned sha256 fbaebb...), DBAASP canonical pull (546 seqs);
    T1 max-Jaccard < 0.35, T2 0.35-0.6, T3 > 0.6 (near-duplicate);
 3. names discovery candidates ONLY from T1/T2 with ensemble p >= 0.7 (top 10);
    T3 hits reported as rediscovery controls;
 4. evidence cards: scores, tiers, physicochemical profile, signal-peptide
    flag + ablated re-score (N-term 25 aa window hydrophobic core heuristic),
    hemolysis-risk flag (heuristic: KD mean > 0.5 AND hydrophobic moment >
    0.6 - documented as a flag, not a validated predictor).
"""
from __future__ import annotations
import json, pathlib
import numpy as np
import torch

from studies.study14_feb2020_full import (rich_features,
    load_fasta_seqs, train_binary, train_binary_shared_adj, MAX_LEN)
from studies.study19_feb2020_v2 import rich_features_v2
from studies.study12_candidates import (net_charge, hydrophobic_moment,
    kmer_jaccard, longest_common_substr, KD, BOMAN)
from pepdesign.models.cnn import PepCNN
from pepdesign.models.gnn import PepGNN
from pepdesign.models.ensemble import predict_scores, predict_scores_shared_adj
from pepdesign.encoding import combined_features
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier

ROOT = pathlib.Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw/feb2020"
WEIGHTS = {"cnn7": .20, "cnn27": .20, "gnn": .15, "rf": .20, "hgb": .25}

def cnn_inputs(seqs, max_len):
    return (torch.stack([torch.from_numpy(combined_features(x, max_len)) for x in seqs]),)

def main():
    amp = load_fasta_seqs(RAW/"AMPS_02182020.fasta")
    dec = load_fasta_seqs(RAW/"DECOYS_02182020.fasta")
    seqs = amp + dec
    y = np.array([1]*len(amp) + [0]*len(dec))
    assert (len(amp), len(dec)) == (2021, 2021)
    print("training final ensemble on full Feb2020 ...", flush=True)
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
        parts = {"cnn7": float(c7), "cnn27": float(c27), "gnn": float(g), "rf": float(r), "hgb": float(h)}
        return sum(WEIGHTS[k]*parts[k] for k in WEIGHTS), parts

    pool = []
    hdrs = {}
    cur_h = None
    for line in open(ROOT/"data/raw/uniprot_uncharacterized_10-60.fasta"):
        line = line.strip()
        if line.startswith(">"):
            cur_h = line[1:]
        elif line and set(line) <= set("ACDEFGHIKLMNPQRSTVWY") and cur_h:
            pool.append(line); hdrs[line] = cur_h
    print(f"pool: {len(pool)} sequences", flush=True)

    refs = {}
    refs["APD6"] = load_fasta_seqs(ROOT/"data/raw/apd6_natural_2024.fasta")
    dr = []
    for i, line in enumerate(open(ROOT/"data/raw/novelty_refs/dramp_natural_amps.txt")):
        if i == 0: continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) > 1 and parts[1] and set(parts[1]) <= set("ACDEFGHIKLMNPQRSTVWY") and 5 <= len(parts[1]) <= 150:
            dr.append(parts[1])
    refs["DRAMP"] = dr
    refs["DBAASP"] = load_fasta_seqs(ROOT/"data/raw/novelty_refs/dbaasp_all.fasta")
    print({k: len(v) for k, v in refs.items()}, flush=True)

    def audit(s):
        best = {"APD6": (0.0, 0), "DRAMP": (0.0, 0), "DBAASP": (0.0, 0)}
        for name, corpus in refs.items():
            jmax, lmax = 0.0, 0
            for r in corpus:
                j = kmer_jaccard(s, r)
                if j > jmax: jmax = j
                if jmax >= 0.35:
                    l = longest_common_substr(s, r)
                    if l > lmax: lmax = l
            best[name] = (round(jmax, 3), lmax)
        overall = max(v[0] for v in best.values())
        tier = "T1" if overall < 0.35 else ("T2" if overall <= 0.6 else "T3")
        return best, overall, tier

    def card(name, s, origin):
        e, parts = ens(s)
        best, overall, tier = audit(s)
        if len(s) >= 10:
            wmeans = [float(np.mean([KD[a] for a in s[i:i+10]])) for i in range(min(20, len(s)-9))]
            sig = max(wmeans) > 1.5
        else:
            sig = False
        abl_e = None
        if sig and len(s) > 30:
            abl_e, _ = ens(s[25:])
        kd_mean = float(np.mean([KD[a] for a in s]))
        hm = hydrophobic_moment(s)
        return {"name": name, "origin": origin, "sequence": s,
                "header": hdrs.get(s),
                "ensemble_p": round(float(e), 4), "branch_scores": {k: round(v,4) for k,v in parts.items()},
                "novelty": {"per_db": {k: {"max_3mer_jaccard": v[0], "lcs": v[1]} for k, v in best.items()},
                            "overall_max_jaccard": round(overall, 3), "tier": tier},
                "physicochemical": {"length": len(s), "net_charge_ph7.4": round(net_charge(s),2),
                                    "kd_mean": round(kd_mean,2), "hydrophobic_moment": round(hm,2),
                                    "boman": round(float(np.mean([BOMAN[a] for a in s])),2)},
                "signal_peptide_flag": bool(sig), "ablated_ensemble_p": (round(float(abl_e),4) if abl_e is not None else None),
                "hemolysis_risk_flag": bool(kd_mean > 0.5 and hm > 0.6),
                "falsifiable_prediction": f"ensemble p(AMP)={e:.3f}; testable by broth-microdilution MIC assay vs E. coli K-12 and S. aureus ATCC 25923"}

    cards, controls = [], []
    scored = []
    for s in pool:
        e, _ = ens(s)
        scored.append((e, s))
    scored.sort(reverse=True)
    print("pool scored; auditing top 40 + study12 candidates", flush=True)
    seen = set()
    for e, s in scored[:40]:
        if s in seen: continue
        seen.add(s)
        c = card(f"pool-{len(cards)+len(controls)+1}", s, "uncharacterized_pool")
        (cards if c["novelty"]["tier"] in ("T1","T2") else controls).append(c)
        print(f"  {c['name']} p={c['ensemble_p']:.3f} tier={c['novelty']['tier']} jac={c['novelty']['overall_max_jaccard']}", flush=True)
    prior = json.load(open(ROOT/"results/study12_candidates.json"))["candidates"]
    prior_cards = []
    for pc in prior:
        c = card(pc["name"], pc["sequence"], "study12_prior_candidate")
        prior_cards.append(c)
        print(f"  {c['name']} p={c['ensemble_p']:.3f} tier={c['novelty']['tier']} jac={c['novelty']['overall_max_jaccard']}", flush=True)

    named = [c for c in cards if c["ensemble_p"] >= 0.7][:10]
    for i, c in enumerate(named, 1):
        c["name"] = f"InstiAMP-21-{i}"
    out = {"protocol": "study21 locked: full-data study19 ensemble re-screen + A2.2 tiered novelty audit",
           "refs": {k: len(v) for k, v in refs.items()},
           "pool_size": len(pool),
           "named_candidates": named,
           "prior_candidates_rescored": prior_cards,
           "rediscovery_controls": controls,
           "top40_tier_counts": {"T1T2": len(cards), "T3": len(controls)}}
    (ROOT/"results/study21_rescreen.json").write_text(json.dumps(out, indent=1))
    print(f"named {len(named)} candidates; controls {len(controls)}; wrote results/study21_rescreen.json", flush=True)

if __name__ == "__main__":
    main()
