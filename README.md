# MEGA27-09: In-Silico Peptide Design Suite

Five benchmarked in-silico peptide-design studies plus a discovery screen, built on real
open datasets (UniProt, IEDB, published benchmark corpora). Core ML: CNN and GNN
architectures (pure-PyTorch graph convolutions), evaluated with cross-validated
AUROC/AUPRC, bootstrap confidence intervals, and permutation tests against published
baselines.

## Studies
1. **AMP-Forge** - antimicrobial peptide design (CNN+GNN ensemble vs AMPlify-style baselines on APD3/UniProt-derived data)
2. **HLA-Bind** - peptide-MHC-I binding prediction (IEDB MHC ligand data)
3. **CPP-Router** - cell-penetrating peptide design (CPPsite 2.0)
4. **SoluPept** - peptide solubility/aggregation prediction (PROSO II corpus)
5. **ACP-Select** - anticancer peptide design (iACP/ACPred corpora)
6. **Discovery screen** - designed-AMP screen against resistance-relevant Gram-negative targets

All CI runs are hermetic (bundled fixtures). Live dataset fetches are marked `live`
and run only outside CI.
