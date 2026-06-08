# sexual-selection-plants

An ovule-gated model of pollen competition and the evolution of post-pollination
prezygotic reproductive isolation (RI) in flowering plants — theory, simulation, and
comparative analysis.

The central hypothesis: where effective pollen load exceeds ovule number
(K = P_eff / O > 1), pollen competition intensifies, maternal pistils evolve stricter
filtering, and prezygotic RI arises as a state-dependent **passenger** of the
reproductive arena — predicted to be **directional** (stronger when the few-ovule,
high-K lineage is the maternal parent) and to **decouple** from the more symmetric,
divergence-tracking postzygotic barriers.

## Layout

```
sexual-selection-plants/
├── data/
│   ├── raw/2025_crossability/   Immutable inputs: Tables 1/2/2b, the Zuntini PI tree,
│   │                            K2P genetic distances, and column descriptions.
│   ├── templates/               Blank input templates (crosses, lineage traits, tree).
│   └── processed/               Derived, regenerable datasets (e.g. the n=110 set).
├── theory/
│   ├── theory_note/             The maths: ovule-gated model + derivation of A_ij.
│   ├── manuscript/              Manuscript skeleton (main.tex) and references.bib.
│   └── notes/                   Working notes and reflective drafts.
├── simulation/
│   ├── src/                     Simulation engine (consolidate to one canonical version).
│   ├── configs/                 Parameter grids (psi, alpha, recombination, pleiotropy).
│   └── results/                 Simulation outputs (gitignored; keep a tiny demo only).
├── analysis/
│   ├── 2025_reproduction/       Reproduction of the team's PGLS/phylolm results.
│   ├── asymmetry/               New A_ij analyses (the interaction & UI-discrimination).
│   └── R/                       R analysis templates (e.g. brms).
├── results/figures/             Publication figures and tables.
├── docs/                        Task lists and project documentation.
└── archive/                     Superseded or unvetted material (kept, not deleted).
```

## Workflow (the agreed sequence)

1. **Theory** (`theory/`) — fix the mechanism before fitting anything. The competition
   ratio K = P_eff/O is the hinge; A_ij is derived in the Sobel–Chen RI index.
2. **Simulation** (`simulation/`) — three layers: (i) ecological competition engine,
   (ii) evolving maternal filters under independent/linked/pleiotropic architectures,
   (iii) realistic data damage on phylogenies, to learn which datasets can carry the claim.
3. **Analysis** (`analysis/`) — apply a method suited to *pairwise, directional* cross
   data (not naive PGLS), with unilateral incompatibility as a named rival from the start.

## Reproducibility notes

- `data/raw/` is treated as read-only. Everything in `data/processed/` is regenerable
  from `data/raw/` by scripts in `analysis/`.
- A single Python environment is pinned in `requirements.txt` at the repo root.
- `archive/ovule_gated_passenger_RI/` is a previously auto-generated package retained
  for reference. Its `results_demo/` was produced wholesale and is **unvetted** — do not
  cite or build on it until each step is re-derived.

## Key caveats carried forward

- The 2025 "PGLS Brownian" numbers match ordinary least squares, not a correctly applied
  phylogenetic regression; RI shows ~zero phylogenetic signal on this tree.
- The ovule–RI signal is confined to mixed-mating pairs and vanishes among outcrossers,
  so analyses should restrict to outcrossers (mating system is a separate lever on K).
- Reciprocal-cross asymmetry and the SI×SC rule (unilateral incompatibility) can be
  collinear; check before interpreting any asymmetry result.
