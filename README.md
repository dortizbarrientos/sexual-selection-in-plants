# sexual-selection-in-plants

An ovule-gated model of pollen competition and the evolution of post-pollination
prezygotic reproductive isolation (RI) in flowering plants — **theory, simulation, and
comparative analysis.**

**Central hypothesis.** Where effective pollen load exceeds ovule number
($K = P_{\text{eff}}/O > 1$), pollen competition intensifies, maternal pistils evolve
stricter filtering, and prezygotic RI arises as a state-dependent **passenger** of the
reproductive arena. The signature is **directional** — RI is stronger when the
few-ovule, high-$K$ lineage is the maternal parent — and **decouples** from the more
symmetric, divergence-tracking postzygotic barriers.

---

## Layout

```
sexual-selection-in-plants/
├── data/
│   ├── raw/2025_crossability/   Immutable inputs: Tables 1/2/2b, Zuntini PI tree,
│   │                            K2P distances, column descriptions. READ-ONLY.
│   ├── templates/               Blank input templates (crosses, lineage traits, tree).
│   └── processed/               Derived, regenerable datasets.
├── theory/
│   ├── theory_note/             The maths: ovule-gated model + derivation of A_ij.
│   ├── manuscript/              Manuscript skeleton (main.tex) + references.bib.
│   └── notes/                   Working notes and reflective drafts.
├── simulation/
│   ├── src/                     Simulation engine (Layer 1 built; Layers 2+ to come).
│   ├── configs/                 Parameter grids (psi, alpha, recombination, pleiotropy).
│   └── results/demo/            Small demo outputs kept in-repo (full runs gitignored).
├── analysis/
<<<<<<< HEAD
│   ├── 00_2025_reproduction/    Reproduction of the team's PGLS/phylolm results (R).
│   ├── 01_2026_analyses/        Data audit + model-readiness assessment (Python).
│   ├── 02_2026_asymmetry/       (empty) home for the new A_ij models.
│   ├── 03_ovule_project_audit/  Outputs of the reciprocal-asymmetry audit.
│   └── 04_reciprocal_family_asymmetry_audit/  Script that produces the 03/ outputs.
=======
│   ├── 00_2025_reproduction/          Reproduction of the team's PGLS/phylolm results (R).
│   ├── 01_2026_data_audit/            Descriptive data audit + model-readiness (Python).
│   ├── 02_reciprocal_asymmetry_audit/ Reciprocal-cross asymmetry audit: script + outputs/.
│   └── 03_asymmetry_models/           (empty) home for the new A_ij models.
>>>>>>> 0324d55 (Tidy analysis/: rename data audit, merge reciprocal-asymmetry script+outputs, add asymmetry-models folder; sync README)
├── results/figures/             Publication figures and tables.
├── docs/                        Task lists and project documentation.
├── scripts/                     One-off maintenance scripts.
└── archive/                     Superseded or unvetted material (kept, not deleted).
```

---

## Code map

Every script, what it does, and whether it is current. Detail in each folder's notes.

| Path | Lang | Purpose | Status |
|---|---|---|---|
| `simulation/src/layer1_competition_engine.py` | Python | **Layer 1 ecological engine.** Draws NegBinomial pollen loads with mean `psi*c*O**alpha`, forms `K = P/O`, the opportunity index `Omega`, and the truncation-selection differential `S_z(K)`. Self-validates `S_z` against direct Monte Carlo and the Falconer intensity table; sweeps the `psi x alpha x O` grid; writes a grid CSV and three figures. No RI yet. | **current, validated** |
| `analysis/00_2025_reproduction/00_PI_ovule_number_analyses_2026_01_10.R` | R | The team's comparative pipeline: assembles ovule numbers from Table 2/2b, merges K2P distances, fits PGLS (`nlme::gls` + `corBrownian`) and `phylolm` (BM/OU/lambda/delta/kappa), tests ovule number (continuous and categorical) and asymmetric RI, and builds the pruned Zuntini tree. | reference |
| `analysis/00_2025_reproduction/R/analysis_brms_template.R` | R | `brms` scaffold for a Bayesian phylogenetic mixed model of RI. Builds a competition proxy `K_proxy` (pollen load if available, else `psi * pollen_production` for a sensitivity `psi`), `logK`/`logO`, and joins maternal/paternal traits onto each cross. Reads the `data_templates/` files; not yet run on real data. | template |
<<<<<<< HEAD
| `analysis/01_2026_analyses/ovule_project_analysis.py` | Python | **Descriptive data audit and model-readiness** (numpy/pandas/matplotlib only — no statistical-model libraries). Coverage, missingness, ovule-number confidence, phylogeny-tip and family coverage, descriptive correlations, RI by measure, K2P-vs-RI, ovule-difference-vs-RI, directional RI by ovule ratio. Emits ~10 figures, ~15 tables (CSV + TeX), `model_ready_PI_pairs.csv`, cleaned datasets, and a LaTeX report. Its RI-vs-predictor panels are exploratory, **not** phylogenetically corrected. | current |
| `analysis/04_reciprocal_family_asymmetry_audit/create_reciprocal_family_audit.py` | Python | Builds the **reciprocal-cross asymmetry audit** that feeds the A_ij analysis (standard library only — `csv`, `collections`): finds species pairs measured in *both* directions, computes the pair-level asymmetry, aggregates to family level, and a "funnel" of how many pairs/families survive the both-directions requirement. Writes the outputs in `analysis/03_ovule_project_audit/`. | current |
=======
| `analysis/01_2026_data_audit/ovule_project_analysis.py` | Python | **Descriptive data audit and model-readiness** (numpy/pandas/matplotlib only — no statistical-model libraries). Coverage, missingness, ovule-number confidence, phylogeny-tip and family coverage, descriptive correlations, RI by measure, K2P-vs-RI, ovule-difference-vs-RI, directional RI by ovule ratio. Emits ~10 figures, ~15 tables (CSV + TeX), `model_ready_PI_pairs.csv`, cleaned datasets, and a LaTeX report. Its RI-vs-predictor panels are exploratory, **not** phylogenetically corrected. | current |
| `analysis/02_reciprocal_asymmetry_audit/create_reciprocal_family_audit.py` | Python | Builds the **reciprocal-cross asymmetry audit** that feeds the A_ij analysis (standard library only — `csv`, `collections`): finds species pairs measured in *both* directions, computes the pair-level asymmetry, aggregates to family level, and a "funnel" of how many pairs/families survive the both-directions requirement. Writes to its own `outputs/` subfolder. | current |
>>>>>>> 0324d55 (Tidy analysis/: rename data audit, merge reciprocal-asymmetry script+outputs, add asymmetry-models folder; sync README)
| `archive/legacy_sims/simulate_ovule_gated_RI.py` | Python | Earlier single-file engine of the full `O -> K -> F -> RI` chain with demo lineages, pairs, and figures. Superseded by the layered rebuild. | archived |
| `archive/legacy_sims/sim1-june-8.{py,R}`, `sim2-june-8.{py,R}` | Py/R | Early simulators that *impose* a linear `RI ~ ovule/distance` signal plus a skew-normal hybrid-vigour tail and fit a Bambi/PyMC phylogenetic mixed model. They assume the pattern rather than generating it; kept only as estimation-side seeds. | archived |
| `archive/ovule_gated_passenger_RI/` | mixed | Original auto-generated package (LaTeX, templates, engine, brms template, `results_demo/`). The demo outputs are **unvetted**; reference only. | archived |
| `scripts/00_reorganize_repo.sh` | shell | One-off repository restructuring helper. | utility |
| `theory/manuscript/main.tex`, `references.bib` | TeX | Manuscript skeleton and bibliography. | skeleton |

Run the engine and point its outputs at the demo folder:

```bash
python simulation/src/layer1_competition_engine.py --outdir simulation/results/demo/layer1
```

---

## Workflow (the agreed sequence)

1. **Theory** (`theory/`) — fix the mechanism before fitting anything. The competition
   ratio $K = P_{\text{eff}}/O$ is the hinge; $A_{ij}$ is derived in the Sobel–Chen RI index.
2. **Simulation** (`simulation/`) — layered, validated one stage at a time:
   (1) ecological competition engine [**done**]; (2) evolving maternal filters under
   independent / linked / pleiotropic architectures; (3) realistic data damage on
   phylogenies, to learn which datasets can carry the claim.
3. **Analysis** (`analysis/`) — a method suited to *pairwise, directional* cross data
   (not naive PGLS), with unilateral incompatibility as a named rival from the start.

---

## Reproducibility notes

- `data/raw/` is read-only; everything in `data/processed/` is regenerable from it by
  scripts in `analysis/`.
- One pinned Python environment in `requirements.txt` at the repo root.
- `archive/` holds superseded or unvetted material, retained for reference. In
  particular `archive/ovule_gated_passenger_RI/results_demo/` was produced wholesale and
  should not be cited or built on until each step is re-derived.

---

## Key caveats carried forward

- The 2025 "PGLS Brownian" numbers match ordinary least squares, not a correctly applied
  phylogenetic regression; RI shows ~zero phylogenetic signal on this tree.
- The ovule–RI signal is confined to mixed-mating pairs and vanishes among outcrossers,
  so analyses should restrict to outcrossers — mating system is a separate lever on $K$.
- Reciprocal-cross asymmetry and the SI×SC rule (unilateral incompatibility) can be
  collinear; check that before interpreting any asymmetry result.
<<<<<<< HEAD
- The `01_2026_analyses` audit is **descriptive only** (no PGLS / phylogenetic model);
=======
- The `01_2026_data_audit` audit is **descriptive only** (no PGLS / phylogenetic model);
>>>>>>> 0324d55 (Tidy analysis/: rename data audit, merge reciprocal-asymmetry script+outputs, add asymmetry-models folder; sync README)
  treat its RI-vs-predictor panels as exploratory and uncorrected, consistent with the
  points above. Its coverage and missingness tables are reliable bookkeeping.

---

## Known issues / portability

<<<<<<< HEAD
- Both 2026 audit scripts (`01_2026_analyses/ovule_project_analysis.py` and
  `04_reciprocal_family_asymmetry_audit/create_reciprocal_family_audit.py`) hard-code
=======
- Both 2026 audit scripts (`01_2026_data_audit/ovule_project_analysis.py` and
  `02_reciprocal_asymmetry_audit/create_reciprocal_family_audit.py`) hard-code
>>>>>>> 0324d55 (Tidy analysis/: rename data audit, merge reciprocal-asymmetry script+outputs, add asymmetry-models folder; sync README)
  their input directory as `/mnt/data/` (the sandbox they were generated in). Repoint
  these to `data/raw/2025_crossability/` before running them in the repo.
- `archive/legacy_sims/simulate_ovule_gated_RI.py` is the most complete legacy
  precursor: it sketches the whole `tree -> ovule states -> K -> filters -> reciprocal
  RI -> data damage` chain in one file. It is superseded by the layered rebuild but is
  the best single reference for Layers 2–6.
- `archive/legacy_sims/sim2-june-8.R` is mislabeled — it contains Python with a stray
  chat preamble line, and a runtime `pip` self-installer. Archived; do not run.
<<<<<<< HEAD

---

=======
>>>>>>> 0324d55 (Tidy analysis/: rename data audit, merge reciprocal-asymmetry script+outputs, add asymmetry-models folder; sync README)
