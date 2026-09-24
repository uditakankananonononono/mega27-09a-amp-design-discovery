# MEGA27 item 9a: antimicrobial peptide design and discovery

Research-grade, reproducible *in silico* AMP classification and candidate
screening with CNN, residue-chain GNN, compositional baselines, 10-mer
cluster-aware splitting, external APD6 validation, a published AMP Scanner
benchmark comparison, and a score-plus-novelty CLI. This repository covers
**item 9 part 1 only**; HLA, CPP, solubility and ACP work are separate projects.

**Evidence boundary:** AMP classifier scores are not calibrated MIC estimates.
No candidate is an experimentally established antibiotic. The strongest
sequence candidate, UniProt P0ACW4/YdcA, has a known N-terminal signal peptide;
removing it drops the frozen ensemble score from 0.908 to 0.630. Synthesis,
activity testing and toxicity testing are needed. On the original Veltri 2018
benchmark our stacked model is near AMP Scanner but below ACEP. Study 13's
updated-dataset CV has a record-filtering flaw and must not be used for a
head-to-head claim; study 14 retains all 4,042 records for an audited
comparison. Published fold-level scores are not available, so a numerical
lead alone does not establish a statistically significant record.

## Reproduce

Create an environment with Python 3.10+, PyTorch, NumPy, scikit-learn,
Biopython, Matplotlib, pytest and install with `pip install -e .`.

```sh
python -m pytest -q tests
python studies/fetch_published_benchmarks.py
python studies/study10_veltri.py
python studies/study11_veltri_stack.py
python studies/study13_feb2020_cv.py  # filtered ablation; not comparable
python studies/study14_feb2020_full.py # all 2,021 AMP + 2,021 decoy
python studies/study12_candidates.py
python studies/study15_signal_ablation.py
python studies/make_sota_section.py
python studies/make_figs2.py
cd paper && pdflatex -interaction=nonstopmode main.tex
```

The large-scale UniProt benchmark and APD6 studies require `fetch_large.py`,
`build_large_dataset.py`, and `fetch_uncharacterized.py` first, then
`study07_amp_scaleup.py`, `study08_discovery.py`, and `study09_apd_validation.py`.
Fetch scripts call live public services; hermetic tests do not. Check
`results/benchmark_sources.json` for the downloaded benchmark FASTA hashes.

## Screen tool

```sh
python -m pepdesign.cli CANDIDATES.fasta \
  --cnn results/study07_cnn.pt --gnn results/study07_gnn.pt \
  --reference data/raw/apd6_natural_2024.fasta --out screen.json
```

The score is an average of CNN and GNN binary classifier outputs. Novelty is
`1 - max Jaccard_3(candidate, reference)` over **individual** reference
sequences, not the fraction of never-seen k-mers or a proof of absence from
all databases. The saved demonstration is `results/tool_screen_demo.json`.

See `paper/main.tex` and `paper/main.pdf` for methods, limitations and result
tables. The PDF is generated after the latest study completes; check its
revision against the corresponding code commit rather than assuming a stale
local copy matches this README.
