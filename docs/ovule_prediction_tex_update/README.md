# Ovule-number reciprocal-asymmetry analysis update

This package contains a TeX note updated with the statistical analysis plan and a reproducible Python script for first-pass analyses.

## Main files

- `ovule_number_pollen_competition_analysis_ready.pdf`: compiled TeX note.
- `ovule_number_pollen_competition_analysis_ready.tex`: updated TeX source.
- `scripts/ovule_asymmetry_analysis.py`: reproducible Python analysis script.
- `data/reciprocal_pair_level_SC.csv`: one row per reciprocal species pair, with exact directional Sobel-Chen RI values and signed ovule contrasts.
- `data/directional_pair_level_SC.csv`: two rows per reciprocal pair, one for each maternal direction.
- `tables/model_summary.csv`: first-pass model coefficients.
- `tables/sign_and_correlation_tests.csv`: raw sign and correlation checks.
- `tables/permutation_test.csv`: within-family ovule-label permutation test.
- `figures/asymmetry_vs_signed_ovule_contrast.png`: raw quadrant diagnostic.

## Sign convention

`A_12 = RI_SC(sp1 mother, sp2 father) - RI_SC(sp2 mother, sp1 father)`

`Delta log O_12 = log(O_2) - log(O_1)`

A positive coefficient means that the fewer-ovule maternal direction has stronger prezygotic reproductive isolation.

## Re-run the analysis

From a directory containing the input CSV and tree files:

```bash
python scripts/ovule_asymmetry_analysis.py \
  --input-dir . \
  --output-dir ovule_prediction_analysis_update \
  --n-permutations 999
```

The script also has an optional flag:

```bash
--use-exact-table2b-fill
```

This uses exact, non-removed, species-level Table2b ovule values to fill missing Table2 ovule numbers. In the current data, the fill candidates do not change the main reciprocal asymmetry-ready sample size.
