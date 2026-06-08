# Ovule-gated pollen competition and passenger reproductive isolation

This folder contains a first working package for the ovule-gated passenger RI model.

## Contents

- `latex/main.tex` — Overleaf-ready theory note.
- `latex/references.bib` — BibTeX references included for later manuscript conversion. The current `main.tex` also has a manual bibliography, so it compiles with pdfLaTeX alone.
- `simulation/simulate_ovule_gated_RI.py` — Python simulation sandbox.
- `simulation/analysis_brms_template.R` — R scaffold for empirical data analysis.
- `data_templates/lineage_traits_template.csv` — Example lineage-level data structure.
- `data_templates/crosses_template.csv` — Example reciprocal-cross data structure.
- `data_templates/phylogeny_template.nwk` — Minimal placeholder Newick tree.
- `results_demo/` — Demo outputs generated from the simulation.

## Compile the LaTeX note

In Overleaf, upload the contents of the `latex` folder and compile `main.tex` with pdfLaTeX. The bibliography is embedded in `main.tex` for immediate compilation; `references.bib` is included as a convenience for later conversion to BibTeX/natbib workflow.

Locally, from the `latex` folder:

```bash
pdflatex main.tex
pdflatex main.tex
```

## Run the simulation

Create a fresh Python environment if you want to keep things clean:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Run a compact demonstration grid:

```bash
python simulation/simulate_ovule_gated_RI.py --outdir results_demo --quick
```

Run a larger grid:

```bash
python simulation/simulate_ovule_gated_RI.py \
  --outdir results_full \
  --n_tips 120 \
  --n_reps 50 \
  --alphas 0.25 0.50 0.75 1.00 1.25 \
  --psis 0.001 0.003 0.01 0.02 0.05 0.10 \
  --architectures independent linked pleiotropic
```

## Main simulation outputs

- `demo_lineages.csv` — simulated species traits: ovules, effective pollen load, K, filter strictness, sweep counts.
- `demo_pairs_full.csv` — simulated reciprocal-cross outcomes for all species pairs.
- `scenario_summary_by_replicate.csv` — recovery diagnostics for each replicate.
- `scenario_summary_aggregated.csv` — mean and standard deviation of recovery diagnostics across replicates.
- `fig_K_vs_ovules.png` — whether K declines with ovule number.
- `fig_asymmetry_vs_delta_logK.png` — whether reciprocal asymmetry tracks differences in K.
- `fig_postRI_vs_divergence.png` — whether postzygotic RI tracks divergence time.

## How to read the key diagnostic

The strongest prediction is directional:

```text
asymmetry = RI(i maternal x j paternal) - RI(j maternal x i paternal)
```

should increase with:

```text
log(K_i) - log(K_j)
```

In the aggregated summary, look for columns containing:

```text
asym_model_beta_delta_logK_mean
```

Positive values mean the model recovered the predicted maternal asymmetry. The ovule-only proxy is:

```text
ovule_proxy_beta_delta_logO_proxy_mean
```

This should be positive when alpha < 1, because lower ovule number implies higher K under sublinear pollen-load scaling.

## First empirical use

Start with the reciprocal-cross table. If raw seed counts and ovule denominators are available, analyse those directly. If only proportions are available, use the proportions but state that denominator information was unavailable. If pollen loads are missing, use sensitivity over `psi` and `alpha`, and treat K as a latent or proxy variable.
