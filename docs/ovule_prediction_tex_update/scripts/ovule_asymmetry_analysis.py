#!/usr/bin/env python3
"""
Ovule-number project: reciprocal asymmetry analysis.

This script builds the analysis-ready datasets for the ovule-number/passenger-RI
model and runs a set of first-pass statistical checks.

The biological prediction is directional:

    A_12 = RI_{1<-2} - RI_{2<-1}
    Delta log O_12 = log(O_2) - log(O_1)

With this sign convention, a positive coefficient means that the species with
fewer ovules shows stronger reproductive isolation when it is the maternal
parent. That is the direct empirical prediction of the maternal-filter model.

Required input files:
    Table1_2026_01_09.csv
    Table2_2025_12_15.csv
    Table2b_2025_11_25.csv
    Genetic_distances_results.csv
    All_PI_pairs_Zuntini_2025_10.tre     optional but recommended

Example:
    python scripts/ovule_asymmetry_analysis.py \
        --input-dir /mnt/data \
        --output-dir /mnt/data/ovule_prediction_analysis_update \
        --n-permutations 1999

The script writes:
    data/reciprocal_pair_level_SC.csv
    data/directional_pair_level_SC.csv
    tables/model_summary.csv
    tables/family_slopes.csv
    tables/leave_one_family_out.csv
    tables/permutation_test.csv
    figures/*.png and *.pdf

The models here are deliberately conservative and transparent. They are not a
replacement for the final Bayesian/phylogenetic model, but they provide the
correct directional data structure and the first stress tests.
"""

from __future__ import annotations

import argparse
import math
import re
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.formula.api as smf

try:
    from Bio import Phylo  # type: ignore
except Exception:  # pragma: no cover - optional dependency
    Phylo = None


# -----------------------------------------------------------------------------
# Small utilities
# -----------------------------------------------------------------------------


def normalise_species_name(x: object) -> str | None:
    """Return a conservative matching key for binomial species names."""
    if pd.isna(x):
        return None
    s = str(x).strip().replace("_", " ")
    s = re.sub(r"\s+", " ", s)
    return s.lower()


def unordered_pair_key(a: object, b: object) -> str:
    """A direction-free key for matching species-pair records."""
    aa = normalise_species_name(a)
    bb = normalise_species_name(b)
    if aa is None or bb is None:
        return ""
    return " || ".join(sorted([aa, bb]))


def as_numeric(series: pd.Series) -> pd.Series:
    """Convert a column to numeric while preserving missing values."""
    return pd.to_numeric(series, errors="coerce")


def one_sided_positive_p(coef: float, se: float) -> float:
    """One-sided p-value for the directional prediction beta > 0."""
    if not np.isfinite(coef) or not np.isfinite(se) or se <= 0:
        return np.nan
    z = coef / se
    return float(1.0 - stats.norm.cdf(z))


def save_current_figure(fig_dir: Path, stem: str) -> None:
    """Save the current matplotlib figure as both PNG and PDF."""
    fig_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(fig_dir / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.savefig(fig_dir / f"{stem}.pdf", bbox_inches="tight")
    plt.close()


def clean_flag_yes_no(x: object, default: str = "no") -> str:
    if pd.isna(x):
        return default
    return str(x).strip().lower()


@dataclass
class Paths:
    input_dir: Path
    output_dir: Path

    @property
    def data_dir(self) -> Path:
        return self.output_dir / "data"

    @property
    def table_dir(self) -> Path:
        return self.output_dir / "tables"

    @property
    def fig_dir(self) -> Path:
        return self.output_dir / "figures"

    @property
    def log_dir(self) -> Path:
        return self.output_dir / "logs"

    def make_dirs(self) -> None:
        for d in [self.output_dir, self.data_dir, self.table_dir, self.fig_dir, self.log_dir]:
            d.mkdir(parents=True, exist_ok=True)


# -----------------------------------------------------------------------------
# Loading and wrangling
# -----------------------------------------------------------------------------


def load_inputs(paths: Paths) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, Optional[object]]:
    """Load all input files. The tree is optional."""
    input_dir = paths.input_dir
    t1 = pd.read_csv(input_dir / "Table1_2026_01_09.csv")
    t2 = pd.read_csv(input_dir / "Table2_2025_12_15.csv")
    t2b = pd.read_csv(input_dir / "Table2b_2025_11_25.csv")
    gd = pd.read_csv(input_dir / "Genetic_distances_results.csv")

    tree = None
    tree_path = input_dir / "All_PI_pairs_Zuntini_2025_10.tre"
    if tree_path.exists() and Phylo is not None:
        try:
            tree = Phylo.read(str(tree_path), "newick")
        except Exception:
            tree = None

    return t1, t2, t2b, gd, tree


def build_ovule_lookup(t2: pd.DataFrame, t2b: pd.DataFrame, table1_species: Iterable[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build the lineage lookup and identify exact Table2b fills.

    Table 2 remains the primary source. If Table 2 has a species row but no
    ovule number, and Table 2b contains an exact, non-removed, species-level
    numeric value for the same binomial, the script records a possible fill.
    The filled value is used in the *_with_exact_table2b_fill columns, but the
    original Table 2 values are kept unchanged.
    """
    t2 = t2.copy()
    t2b = t2b.copy()
    t2["species_key"] = t2["species"].map(normalise_species_name)
    t2b["species_key"] = t2b["species"].map(normalise_species_name)
    t2["ovule_number"] = as_numeric(t2["ovule_number"])
    t2["confidence"] = as_numeric(t2["confidence"])
    t2b["number"] = as_numeric(t2b["number"])
    t2b["remove_flag"] = t2b["Remove"].map(lambda x: clean_flag_yes_no(x)).eq("yes")

    lookup = t2[["genus", "species", "species_key", "ovule_number", "confidence", "life_history", "mating_system"]].copy()
    lookup["ovule_number_filled"] = lookup["ovule_number"]
    lookup["ovule_source"] = np.where(lookup["ovule_number"].notna(), "Table2", "missing")
    lookup["filled_from_table2b"] = False

    # Candidate exact species-level values from Table2b.
    candidates = t2b.loc[
        (~t2b["remove_flag"]) &
        (t2b["number"].notna()) &
        (t2b["level"].astype(str).str.lower().str.strip().eq("species")),
        ["family", "genus", "species", "species_key", "level", "number", "Remove"]
    ].copy()
    candidates = candidates.sort_values(["species_key", "number"])
    exact_best = candidates.groupby("species_key", as_index=False).agg(
        table2b_candidate_count=("number", "size"),
        table2b_number_mean=("number", "mean"),
        table2b_number_min=("number", "min"),
        table2b_number_max=("number", "max"),
    )

    table1_keys = {normalise_species_name(x) for x in table1_species if pd.notna(x)}
    fill_candidates = lookup.merge(exact_best, on="species_key", how="left")
    fill_candidates = fill_candidates.loc[
        fill_candidates["species_key"].isin(table1_keys) &
        fill_candidates["ovule_number"].isna() &
        fill_candidates["table2b_number_mean"].notna() &
        np.isclose(fill_candidates["table2b_number_min"], fill_candidates["table2b_number_max"], equal_nan=False)
    ].copy()

    if not fill_candidates.empty:
        key_to_fill = fill_candidates.set_index("species_key")["table2b_number_mean"].to_dict()
        mask = lookup["species_key"].isin(key_to_fill.keys()) & lookup["ovule_number_filled"].isna()
        lookup.loc[mask, "ovule_number_filled"] = lookup.loc[mask, "species_key"].map(key_to_fill)
        lookup.loc[mask, "ovule_source"] = "Table2b_exact_species_fill"
        lookup.loc[mask, "filled_from_table2b"] = True

    return lookup, fill_candidates


def collapse_genetic_distances(gd: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Collapse K2P distances to one value per unordered species pair."""
    gd = gd.copy()
    gd["species1_key"] = gd["Species_1"].map(normalise_species_name)
    gd["species2_key"] = gd["Species_2"].map(normalise_species_name)
    gd["pair_key"] = [unordered_pair_key(a, b) for a, b in zip(gd["Species_1"], gd["Species_2"])]
    gd["K2P_Distance"] = as_numeric(gd["K2P_Distance"])

    collapsed = gd.groupby("pair_key", as_index=False).agg(
        gd_rows=("K2P_Distance", "size"),
        gd_nonmissing=("K2P_Distance", lambda x: int(x.notna().sum())),
        K2P_Distance=("K2P_Distance", "mean"),
        K2P_min=("K2P_Distance", "min"),
        K2P_max=("K2P_Distance", "max"),
    )
    collapsed["gd_duplicate_pair"] = collapsed["gd_rows"] > 1
    collapsed["gd_range"] = collapsed["K2P_max"] - collapsed["K2P_min"]

    audit = pd.DataFrame({
        "metric": [
            "genetic_distance_rows",
            "unique_unordered_pairs",
            "duplicate_unordered_pairs",
            "rows_with_nonmissing_K2P",
            "rows_with_missing_K2P",
        ],
        "value": [
            len(gd),
            collapsed["pair_key"].nunique(),
            int(collapsed["gd_duplicate_pair"].sum()),
            int(gd["K2P_Distance"].notna().sum()),
            int(gd["K2P_Distance"].isna().sum()),
        ],
    })
    return collapsed, audit


def make_pair_level_table(
    t1: pd.DataFrame,
    ovule_lookup: pd.DataFrame,
    gd_collapsed: pd.DataFrame,
    use_filled_ovules: bool = False,
) -> pd.DataFrame:
    """Create one analysis row per retained species pair."""
    df = t1.copy()
    df["row_id_original"] = np.arange(len(df))
    df["remove_clean"] = df["remove"].map(lambda x: clean_flag_yes_no(x))
    df["retained"] = df["remove_clean"].eq("no")
    df["is_PI_from_PI_name"] = df["PI_name"].notna()
    df["species1_key"] = df["species1"].map(normalise_species_name)
    df["species2_key"] = df["species2"].map(normalise_species_name)
    df["pair_key"] = [unordered_pair_key(a, b) for a, b in zip(df["species1"], df["species2"])]

    for c in [
        "intraspecific_value_sp1", "intraspecific_value_sp2",
        "interspecific_value_sp1_sp2", "interspecific_value_sp2_sp1",
        "sp1_hybrid_proportion", "sp2_hybrid_proportion", "RI",
    ]:
        if c in df.columns:
            df[c] = as_numeric(df[c])

    lookup = ovule_lookup.set_index("species_key")
    ovule_column = "ovule_number_filled" if use_filled_ovules else "ovule_number"
    suffix_source = "filled" if use_filled_ovules else "table2"

    for side in [1, 2]:
        key = f"species{side}_key"
        df[f"ovule_sp{side}"] = df[key].map(lookup[ovule_column])
        df[f"ovule_source_sp{side}"] = df[key].map(lookup["ovule_source"])
        df[f"confidence_sp{side}"] = df[key].map(lookup["confidence"])
        df[f"life_history_sp{side}"] = df[key].map(lookup["life_history"])
        df[f"mating_system_sp{side}"] = df[key].map(lookup["mating_system"])
        df[f"ovule_lookup_source_sp{side}"] = suffix_source

    df = df.merge(gd_collapsed, on="pair_key", how="left")

    # Directional Sobel-Chen indices.
    C11 = df["intraspecific_value_sp1"]
    C22 = df["intraspecific_value_sp2"]
    H12 = df["interspecific_value_sp1_sp2"]
    H21 = df["interspecific_value_sp2_sp1"]
    den12 = C11 + H12
    den21 = C22 + H21

    df["sobel_chen_RI_sp1_mother"] = np.where(den12 > 0, (C11 - H12) / den12, np.nan)
    df["sobel_chen_RI_sp2_mother"] = np.where(den21 > 0, (C22 - H21) / den21, np.nan)
    df["sobel_chen_asymmetry_12"] = df["sobel_chen_RI_sp1_mother"] - df["sobel_chen_RI_sp2_mother"]

    # Older 1 - H/C diagnostic, kept for comparison only.
    df["diagnostic_RI_sp1_mother_ratio"] = np.where(C11 > 0, 1 - H12 / C11, np.nan)
    df["diagnostic_RI_sp2_mother_ratio"] = np.where(C22 > 0, 1 - H21 / C22, np.nan)
    df["diagnostic_asymmetry_ratio_index"] = df["diagnostic_RI_sp1_mother_ratio"] - df["diagnostic_RI_sp2_mother_ratio"]

    # Ovule contrasts. Positive delta means species 1 has fewer ovules than species 2.
    valid_ovules = (df["ovule_sp1"] > 0) & (df["ovule_sp2"] > 0)
    df["log_ovules_sp1"] = np.where(valid_ovules, np.log(df["ovule_sp1"]), np.nan)
    df["log_ovules_sp2"] = np.where(valid_ovules, np.log(df["ovule_sp2"]), np.nan)
    df["log10_ovules_sp1"] = np.where(valid_ovules, np.log10(df["ovule_sp1"]), np.nan)
    df["log10_ovules_sp2"] = np.where(valid_ovules, np.log10(df["ovule_sp2"]), np.nan)
    df["delta_log_ovules_12"] = df["log_ovules_sp2"] - df["log_ovules_sp1"]
    df["delta_log10_ovules_12"] = df["log10_ovules_sp2"] - df["log10_ovules_sp1"]
    df["abs_delta_log_ovules"] = np.abs(df["delta_log_ovules_12"])
    df["nonzero_ovule_difference"] = df["delta_log_ovules_12"].abs() > 1e-12

    # The divergence-scaled prediction.
    df["K2P_delta_log_ovules_12"] = df["K2P_Distance"] * df["delta_log_ovules_12"]

    # Mating-system codings. These are rough modifiers, not SI status.
    outcross_score = {"selfing": 0.0, "mixed-mating": 0.5, "outcrossing": 1.0}
    df["outcross_score_sp1"] = df["mating_system_sp1"].map(outcross_score)
    df["outcross_score_sp2"] = df["mating_system_sp2"].map(outcross_score)
    df["mean_outcross_score"] = (df["outcross_score_sp1"] + df["outcross_score_sp2"]) / 2
    df["delta_outcross_score_12"] = df["outcross_score_sp1"] - df["outcross_score_sp2"]
    df["same_mating_system"] = df["mating_system_sp1"].eq(df["mating_system_sp2"])

    # Readiness flags.
    df["has_reciprocal_interspecific"] = H12.notna() & H21.notna()
    df["has_both_intraspecific_controls"] = C11.notna() & C22.notna()
    df["has_directional_SC_RI"] = df["sobel_chen_RI_sp1_mother"].notna() & df["sobel_chen_RI_sp2_mother"].notna()
    df["has_both_ovules"] = valid_ovules
    df["asymmetry_model_ready"] = (
        df["retained"] & df["has_reciprocal_interspecific"] &
        df["has_both_intraspecific_controls"] & df["has_directional_SC_RI"] &
        df["has_both_ovules"]
    )
    df["both_confidence_le_2"] = (
        df["confidence_sp1"].notna() & df["confidence_sp2"].notna() &
        df["confidence_sp1"].le(2) & df["confidence_sp2"].le(2)
    )
    df["strict_PI_K2P_asymmetry_ready"] = (
        df["asymmetry_model_ready"] & df["is_PI_from_PI_name"] & df["K2P_Distance"].notna()
    )
    df["prediction_sign_agreement"] = (
        np.sign(df["sobel_chen_asymmetry_12"]) == np.sign(df["delta_log_ovules_12"])
    ) & df["nonzero_ovule_difference"] & (df["sobel_chen_asymmetry_12"].abs() > 1e-12)

    return df


def make_directional_table(pair_df: pd.DataFrame) -> pd.DataFrame:
    """Convert pair rows to two directional rows for pair-fixed-effect models."""
    rows: list[dict[str, object]] = []
    ready = pair_df.loc[pair_df["asymmetry_model_ready"]].copy()
    for _, r in ready.iterrows():
        rows.append({
            "pair_key": r["pair_key"],
            "family": r["family"],
            "genus": r["genus"],
            "measure": r["measure"],
            "maternal_direction": "species1_mother",
            "maternal_species": r["species1"],
            "paternal_species": r["species2"],
            "maternal_species_key": r["species1_key"],
            "paternal_species_key": r["species2_key"],
            "RI_SC": r["sobel_chen_RI_sp1_mother"],
            "maternal_ovules": r["ovule_sp1"],
            "paternal_ovules": r["ovule_sp2"],
            "maternal_log_ovules": r["log_ovules_sp1"],
            "negative_maternal_log_ovules": -r["log_ovules_sp1"],
            "maternal_confidence": r["confidence_sp1"],
            "paternal_confidence": r["confidence_sp2"],
            "maternal_mating_system": r["mating_system_sp1"],
            "paternal_mating_system": r["mating_system_sp2"],
            "K2P_Distance": r["K2P_Distance"],
            "is_PI_from_PI_name": r["is_PI_from_PI_name"],
        })
        rows.append({
            "pair_key": r["pair_key"],
            "family": r["family"],
            "genus": r["genus"],
            "measure": r["measure"],
            "maternal_direction": "species2_mother",
            "maternal_species": r["species2"],
            "paternal_species": r["species1"],
            "maternal_species_key": r["species2_key"],
            "paternal_species_key": r["species1_key"],
            "RI_SC": r["sobel_chen_RI_sp2_mother"],
            "maternal_ovules": r["ovule_sp2"],
            "paternal_ovules": r["ovule_sp1"],
            "maternal_log_ovules": r["log_ovules_sp2"],
            "negative_maternal_log_ovules": -r["log_ovules_sp2"],
            "maternal_confidence": r["confidence_sp2"],
            "paternal_confidence": r["confidence_sp1"],
            "maternal_mating_system": r["mating_system_sp2"],
            "paternal_mating_system": r["mating_system_sp1"],
            "K2P_Distance": r["K2P_Distance"],
            "is_PI_from_PI_name": r["is_PI_from_PI_name"],
        })
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Summaries
# -----------------------------------------------------------------------------


def make_funnel_table(pair_df: pd.DataFrame) -> pd.DataFrame:
    """Data attrition from raw rows to model-ready subsets."""
    conditions = [
        ("All Table 1 rows", pair_df.index.notna()),
        ("Retained rows: remove == no", pair_df["retained"]),
        ("Both reciprocal interspecific directions", pair_df["retained"] & pair_df["has_reciprocal_interspecific"]),
        ("Directional Sobel-Chen RI possible", pair_df["retained"] & pair_df["has_reciprocal_interspecific"] & pair_df["has_both_intraspecific_controls"] & pair_df["has_directional_SC_RI"]),
        ("Asymmetry model-ready: plus both ovule numbers", pair_df["asymmetry_model_ready"]),
        ("Asymmetry-ready with nonzero ovule contrast", pair_df["asymmetry_model_ready"] & pair_df["nonzero_ovule_difference"]),
        ("High-confidence ovules: both confidence <= 2", pair_df["asymmetry_model_ready"] & pair_df["both_confidence_le_2"]),
        ("Strict PI/K2P asymmetry subset", pair_df["strict_PI_K2P_asymmetry_ready"]),
        ("Strict PI/K2P plus confidence <= 2", pair_df["strict_PI_K2P_asymmetry_ready"] & pair_df["both_confidence_le_2"]),
    ]
    rows = []
    for label, mask in conditions:
        m = pd.Series(mask, index=pair_df.index).fillna(False)
        sub = pair_df.loc[m]
        rows.append({
            "condition": label,
            "pair_rows": int(len(sub)),
            "families": int(sub["family"].nunique()) if "family" in sub else 0,
            "species": int(pd.concat([sub["species1"], sub["species2"]]).nunique()) if len(sub) else 0,
        })
    return pd.DataFrame(rows)


def make_family_summary(pair_df: pd.DataFrame) -> pd.DataFrame:
    """Family-level readiness counts."""
    grouped = pair_df.loc[pair_df["retained"]].groupby("family", dropna=False)
    rows = []
    for fam, g in grouped:
        ready = g["asymmetry_model_ready"]
        high = ready & g["both_confidence_le_2"]
        strict = g["strict_PI_K2P_asymmetry_ready"]
        nonzero = ready & g["nonzero_ovule_difference"]
        rows.append({
            "family": fam,
            "retained_rows": int(len(g)),
            "reciprocal_interspecific_pairs": int((g["has_reciprocal_interspecific"]).sum()),
            "directional_SC_RI_pairs": int((g["has_reciprocal_interspecific"] & g["has_both_intraspecific_controls"] & g["has_directional_SC_RI"]).sum()),
            "asymmetry_ready_pairs": int(ready.sum()),
            "asymmetry_ready_nonzero_delta_logO": int(nonzero.sum()),
            "asymmetry_ready_high_confidence": int(high.sum()),
            "strict_PI_K2P_asymmetry_ready": int(strict.sum()),
            "distinct_delta_logO_values_ready": int(g.loc[ready, "delta_log_ovules_12"].round(12).nunique()),
            "species_in_ready_pairs": int(pd.concat([g.loc[ready, "species1"], g.loc[ready, "species2"]]).nunique()) if ready.any() else 0,
            "genera_in_ready_pairs": int(g.loc[ready, "genus"].nunique()) if ready.any() else 0,
        })
    out = pd.DataFrame(rows)
    out = out.sort_values(["asymmetry_ready_pairs", "reciprocal_interspecific_pairs"], ascending=False)
    return out


def make_family_modelability_table(family_summary: pd.DataFrame) -> pd.DataFrame:
    """How many families can support within-family slopes under simple thresholds."""
    rows = []
    for min_pairs in [2, 3, 5, 10]:
        all_mask = (family_summary["asymmetry_ready_pairs"] >= min_pairs) & (family_summary["distinct_delta_logO_values_ready"] >= 2)
        high_mask = (family_summary["asymmetry_ready_high_confidence"] >= min_pairs) & (family_summary["distinct_delta_logO_values_ready"] >= 2)
        rows.append({
            "minimum_pairs_per_family": min_pairs,
            "families_all_numeric_ovules": int(all_mask.sum()),
            "families_high_confidence_proxy": int(high_mask.sum()),
        })
    return pd.DataFrame(rows)


def audit_tree(tree: Optional[object], pair_df: pd.DataFrame) -> pd.DataFrame:
    """Simple tree/tip audit. Does not build a covariance matrix."""
    if tree is None:
        return pd.DataFrame({"metric": ["tree_loaded"], "value": ["no"]})
    tips = [tip.name for tip in tree.get_terminals() if tip.name is not None]
    tip_counts = pd.Series(tips).value_counts()
    pi_names = set(pair_df.loc[pair_df["retained"] & pair_df["is_PI_from_PI_name"], "PI_name"].dropna())
    tip_set = set(tips)
    duplicated = sorted(tip_counts.loc[tip_counts > 1].index.tolist())
    return pd.DataFrame({
        "metric": [
            "tree_loaded",
            "tip_instances",
            "unique_tip_labels",
            "duplicated_tip_labels",
            "retained_PI_labels",
            "PI_labels_found_in_tree",
            "PI_labels_missing_from_tree",
        ],
        "value": [
            "yes",
            len(tips),
            len(tip_set),
            "; ".join(duplicated) if duplicated else "none",
            len(pi_names),
            len(pi_names & tip_set),
            len(pi_names - tip_set),
        ],
    })


# -----------------------------------------------------------------------------
# Models
# -----------------------------------------------------------------------------


def fit_ols_model(name: str, data: pd.DataFrame, formula: str, target: str, cov_type: str = "HC3") -> dict[str, object]:
    """Fit an OLS model and return the target coefficient."""
    try:
        fit = smf.ols(formula, data=data).fit(cov_type=cov_type)
        coef = float(fit.params.get(target, np.nan))
        se = float(fit.bse.get(target, np.nan))
        p_two = float(fit.pvalues.get(target, np.nan))
        p_one = one_sided_positive_p(coef, se)
        return {
            "model": name,
            "engine": f"statsmodels OLS, cov_type={cov_type}",
            "n": int(fit.nobs),
            "target_term": target,
            "coefficient": coef,
            "standard_error": se,
            "p_two_sided": p_two,
            "p_one_sided_positive": p_one,
            "r_squared": float(getattr(fit, "rsquared", np.nan)),
            "formula": formula,
            "status": "ok",
        }
    except Exception as e:
        return {
            "model": name,
            "engine": f"statsmodels OLS, cov_type={cov_type}",
            "n": int(len(data)),
            "target_term": target,
            "coefficient": np.nan,
            "standard_error": np.nan,
            "p_two_sided": np.nan,
            "p_one_sided_positive": np.nan,
            "r_squared": np.nan,
            "formula": formula,
            "status": f"failed: {type(e).__name__}: {e}",
        }


def fit_directional_pair_fixed_effect(long_df: pd.DataFrame) -> dict[str, object]:
    """Fit RI(direction) ~ -log maternal ovules + pair fixed effects.

    This is the same prediction as the asymmetry regression, written in a
    long format. Pair fixed effects absorb all symmetric differences between
    species pairs. Standard errors are clustered by pair.
    """
    formula = "RI_SC ~ negative_maternal_log_ovules + C(pair_key)"
    target = "negative_maternal_log_ovules"
    try:
        fit = smf.ols(formula, data=long_df).fit(
            cov_type="cluster",
            cov_kwds={"groups": long_df["pair_key"]},
        )
        coef = float(fit.params.get(target, np.nan))
        se = float(fit.bse.get(target, np.nan))
        return {
            "model": "directional_pair_fixed_effect",
            "engine": "statsmodels OLS with pair fixed effects; SE clustered by pair",
            "n": int(fit.nobs),
            "target_term": target,
            "coefficient": coef,
            "standard_error": se,
            "p_two_sided": float(fit.pvalues.get(target, np.nan)),
            "p_one_sided_positive": one_sided_positive_p(coef, se),
            "r_squared": float(getattr(fit, "rsquared", np.nan)),
            "formula": formula,
            "status": "ok",
        }
    except Exception as e:
        return {
            "model": "directional_pair_fixed_effect",
            "engine": "statsmodels OLS with pair fixed effects; SE clustered by pair",
            "n": int(len(long_df)),
            "target_term": target,
            "coefficient": np.nan,
            "standard_error": np.nan,
            "p_two_sided": np.nan,
            "p_one_sided_positive": np.nan,
            "r_squared": np.nan,
            "formula": formula,
            "status": f"failed: {type(e).__name__}: {e}",
        }




def fit_pair_fixed_effect_fast(pair_df: pd.DataFrame) -> dict[str, object]:
    """Fast equivalent of the two-row pair-fixed-effect model.

    In the long model, RI_SC is regressed on -log maternal ovules with a
    fixed effect for each species pair. Demeaning the two directions within
    each pair reduces the model to a no-intercept regression of the directional
    asymmetry A_12 on Delta log O_12. This avoids hundreds of dummy variables
    and gives the same slope for the ovule term.
    """
    data = pair_df.loc[pair_df["asymmetry_model_ready"] & pair_df["nonzero_ovule_difference"]].copy()
    formula = "sobel_chen_asymmetry_12 ~ 0 + delta_log_ovules_12"
    target = "delta_log_ovules_12"
    try:
        fit = smf.ols(formula, data=data).fit(cov_type="HC3")
        coef = float(fit.params.get(target, np.nan))
        se = float(fit.bse.get(target, np.nan))
        return {
            "model": "directional_pair_fixed_effect_fast",
            "engine": "fast paired-difference equivalent of RI_SC ~ -log maternal ovules + pair fixed effects",
            "n": int(fit.nobs),
            "target_term": target,
            "coefficient": coef,
            "standard_error": se,
            "p_two_sided": float(fit.pvalues.get(target, np.nan)),
            "p_one_sided_positive": one_sided_positive_p(coef, se),
            "r_squared": float(getattr(fit, "rsquared", np.nan)),
            "formula": formula,
            "status": "ok",
        }
    except Exception as e:
        return {
            "model": "directional_pair_fixed_effect_fast",
            "engine": "fast paired-difference equivalent of RI_SC ~ -log maternal ovules + pair fixed effects",
            "n": int(len(data)),
            "target_term": target,
            "coefficient": np.nan,
            "standard_error": np.nan,
            "p_two_sided": np.nan,
            "p_one_sided_positive": np.nan,
            "r_squared": np.nan,
            "formula": formula,
            "status": f"failed: {type(e).__name__}: {e}",
        }


def fit_family_mixed_model(data: pd.DataFrame) -> dict[str, object]:
    """Fit a random family intercept/slope model if it converges.

    Statsmodels MixedLM is a useful diagnostic here, but small families and
    unbalanced data can cause boundary fits. Treat this as a first-pass check,
    not as the final inferential model.
    """
    formula = "sobel_chen_asymmetry_12 ~ delta_log_ovules_12 + C(measure)"
    target = "delta_log_ovules_12"
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            fit = smf.mixedlm(
                formula,
                data=data,
                groups=data["family"],
                re_formula="~delta_log_ovules_12",
            ).fit(method="lbfgs", reml=False, maxiter=1000, disp=False)
        coef = float(fit.params.get(target, np.nan))
        se = float(fit.bse.get(target, np.nan))
        warning_text = " | ".join(sorted({str(w.message) for w in caught}))
        status = "ok" if getattr(fit, "converged", False) else "not_converged"
        if warning_text:
            status = status + "; warnings: " + warning_text[:500]
        return {
            "model": "family_random_intercept_slope",
            "engine": "statsmodels MixedLM, family random intercept and slope",
            "n": int(fit.nobs),
            "target_term": target,
            "coefficient": coef,
            "standard_error": se,
            "p_two_sided": float(fit.pvalues.get(target, np.nan)),
            "p_one_sided_positive": one_sided_positive_p(coef, se),
            "r_squared": np.nan,
            "formula": formula + " + (1 + delta_log_ovules_12 | family)",
            "status": status,
        }
    except Exception as e:
        return {
            "model": "family_random_intercept_slope",
            "engine": "statsmodels MixedLM, family random intercept and slope",
            "n": int(len(data)),
            "target_term": target,
            "coefficient": np.nan,
            "standard_error": np.nan,
            "p_two_sided": np.nan,
            "p_one_sided_positive": np.nan,
            "r_squared": np.nan,
            "formula": formula + " + (1 + delta_log_ovules_12 | family)",
            "status": f"failed: {type(e).__name__}: {e}",
        }


def run_model_suite(pair_df: pd.DataFrame, directional_df: pd.DataFrame) -> pd.DataFrame:
    """Fit first-pass models for the prediction."""
    d_all = pair_df.loc[pair_df["asymmetry_model_ready"]].copy()
    d_high = d_all.loc[d_all["both_confidence_le_2"]].copy()
    d_strict = d_all.loc[d_all["strict_PI_K2P_asymmetry_ready"]].copy()

    rows: list[dict[str, object]] = []
    rows.append(fit_ols_model(
        "broad_raw_slope",
        d_all,
        "sobel_chen_asymmetry_12 ~ delta_log_ovules_12",
        "delta_log_ovules_12",
    ))
    rows.append(fit_ols_model(
        "broad_measure_adjusted",
        d_all,
        "sobel_chen_asymmetry_12 ~ delta_log_ovules_12 + C(measure)",
        "delta_log_ovules_12",
    ))
    rows.append(fit_ols_model(
        "broad_family_and_measure_fixed_effects",
        d_all,
        "sobel_chen_asymmetry_12 ~ delta_log_ovules_12 + C(measure) + C(family)",
        "delta_log_ovules_12",
        cov_type="HC1",
    ))
    rows.append(fit_ols_model(
        "high_confidence_family_and_measure_fixed_effects",
        d_high,
        "sobel_chen_asymmetry_12 ~ delta_log_ovules_12 + C(measure) + C(family)",
        "delta_log_ovules_12",
        cov_type="HC1",
    ))
    rows.append(fit_ols_model(
        "strict_PI_K2P_divergence_interaction",
        d_strict,
        "sobel_chen_asymmetry_12 ~ delta_log_ovules_12 + K2P_Distance + K2P_delta_log_ovules_12 + C(measure)",
        "K2P_delta_log_ovules_12",
    ))
    rows.append(fit_ols_model(
        "mating_system_modifier",
        d_all.loc[d_all["mean_outcross_score"].notna()].copy(),
        "sobel_chen_asymmetry_12 ~ delta_log_ovules_12 * mean_outcross_score + C(measure)",
        "delta_log_ovules_12:mean_outcross_score",
    ))

    for measure, g in d_all.groupby("measure"):
        if len(g) >= 10 and g["delta_log_ovules_12"].nunique() >= 2:
            rows.append(fit_ols_model(
                f"measure_stratum_{measure}",
                g,
                "sobel_chen_asymmetry_12 ~ delta_log_ovules_12",
                "delta_log_ovules_12",
            ))

    rows.append(fit_pair_fixed_effect_fast(pair_df))
    # The final paper should use a partial-pooling family model.
    # Statsmodels MixedLM is fragile for these highly unbalanced family data,
    # so this executable first-pass suite does not run it by default. The
    # fit_family_mixed_model() function above is kept as a diagnostic recipe.

    return pd.DataFrame(rows)


def make_family_slopes(pair_df: pd.DataFrame) -> pd.DataFrame:
    """Fit simple within-family slopes for descriptive plotting."""
    d = pair_df.loc[pair_df["asymmetry_model_ready"]].copy()
    rows = []
    for family, g in d.groupby("family"):
        n_delta = g["delta_log_ovules_12"].round(12).nunique()
        if len(g) >= 3 and n_delta >= 2:
            try:
                fit = smf.ols("sobel_chen_asymmetry_12 ~ delta_log_ovules_12", data=g).fit()
                coef = float(fit.params["delta_log_ovules_12"])
                se = float(fit.bse["delta_log_ovules_12"])
                rows.append({
                    "family": family,
                    "n_pairs": int(len(g)),
                    "distinct_delta_logO": int(n_delta),
                    "slope": coef,
                    "standard_error": se,
                    "p_two_sided": float(fit.pvalues["delta_log_ovules_12"]),
                    "p_one_sided_positive": one_sided_positive_p(coef, se),
                    "mean_asymmetry": float(g["sobel_chen_asymmetry_12"].mean()),
                    "mean_delta_logO": float(g["delta_log_ovules_12"].mean()),
                })
            except Exception as e:
                rows.append({
                    "family": family,
                    "n_pairs": int(len(g)),
                    "distinct_delta_logO": int(n_delta),
                    "slope": np.nan,
                    "standard_error": np.nan,
                    "p_two_sided": np.nan,
                    "p_one_sided_positive": np.nan,
                    "mean_asymmetry": float(g["sobel_chen_asymmetry_12"].mean()),
                    "mean_delta_logO": float(g["delta_log_ovules_12"].mean()),
                    "status": f"failed: {type(e).__name__}: {e}",
                })
    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values("n_pairs", ascending=False)
    return out


def make_leave_one_family_out(pair_df: pd.DataFrame) -> pd.DataFrame:
    """Re-fit the broad measure-adjusted model after dropping one family at a time."""
    d = pair_df.loc[pair_df["asymmetry_model_ready"]].copy()
    rows = []
    for family in sorted(d["family"].dropna().unique()):
        sub = d.loc[d["family"] != family].copy()
        result = fit_ols_model(
            f"drop_{family}",
            sub,
            "sobel_chen_asymmetry_12 ~ delta_log_ovules_12 + C(measure)",
            "delta_log_ovules_12",
        )
        result["dropped_family"] = family
        result["n_pairs_dropped"] = int((d["family"] == family).sum())
        rows.append(result)
    return pd.DataFrame(rows)


def make_sign_tests(pair_df: pd.DataFrame) -> pd.DataFrame:
    """Simple sign and correlation diagnostics."""
    rows = []
    subsets = {
        "all_asymmetry_ready": pair_df["asymmetry_model_ready"],
        "high_confidence": pair_df["asymmetry_model_ready"] & pair_df["both_confidence_le_2"],
        "strict_PI_K2P": pair_df["strict_PI_K2P_asymmetry_ready"],
        "nonzero_ovule_difference": pair_df["asymmetry_model_ready"] & pair_df["nonzero_ovule_difference"],
    }
    for name, mask in subsets.items():
        g = pair_df.loc[mask].copy()
        nonzero = g.loc[(g["sobel_chen_asymmetry_12"].abs() > 1e-12) & (g["delta_log_ovules_12"].abs() > 1e-12)].copy()
        if len(nonzero) > 0:
            successes = int((np.sign(nonzero["sobel_chen_asymmetry_12"]) == np.sign(nonzero["delta_log_ovules_12"])).sum())
            sign_n = int(len(nonzero))
            binom_p = float(stats.binomtest(successes, sign_n, p=0.5, alternative="greater").pvalue)
        else:
            successes, sign_n, binom_p = 0, 0, np.nan
        if len(g) >= 3 and g["delta_log_ovules_12"].nunique() >= 2:
            pearson = stats.pearsonr(g["delta_log_ovules_12"], g["sobel_chen_asymmetry_12"])
            spearman = stats.spearmanr(g["delta_log_ovules_12"], g["sobel_chen_asymmetry_12"], nan_policy="omit")
            pearson_r, pearson_p = float(pearson.statistic), float(pearson.pvalue)
            spearman_r, spearman_p = float(spearman.statistic), float(spearman.pvalue)
        else:
            pearson_r = pearson_p = spearman_r = spearman_p = np.nan
        rows.append({
            "subset": name,
            "n_pairs": int(len(g)),
            "n_sign_test": sign_n,
            "sign_agreements": successes,
            "sign_agreement_fraction": successes / sign_n if sign_n else np.nan,
            "binomial_p_one_sided_positive": binom_p,
            "pearson_r": pearson_r,
            "pearson_p_two_sided": pearson_p,
            "spearman_rho": spearman_r,
            "spearman_p_two_sided": spearman_p,
            "mean_asymmetry": float(g["sobel_chen_asymmetry_12"].mean()) if len(g) else np.nan,
            "sd_asymmetry": float(g["sobel_chen_asymmetry_12"].std()) if len(g) else np.nan,
        })
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Permutation test
# -----------------------------------------------------------------------------


def residualize_against_controls(x: np.ndarray, q_matrix: np.ndarray) -> np.ndarray:
    """Residualise a vector against a control design with orthonormal Q columns."""
    return x - q_matrix @ (q_matrix.T @ x)


def build_control_q(data: pd.DataFrame) -> np.ndarray:
    """Build an orthonormal basis for intercept + family + measure controls."""
    controls = pd.get_dummies(data[["family", "measure"]], drop_first=True, dtype=float)
    X = np.column_stack([np.ones(len(data)), controls.to_numpy(dtype=float)])
    q, _ = np.linalg.qr(X)
    return q


def family_stratified_ovule_permutation(
    pair_df: pd.DataFrame,
    n_permutations: int,
    seed: int = 20260609,
) -> pd.DataFrame:
    """Shuffle ovule values among species within families and re-fit the slope.

    The cross network, families, measure types, and RI values stay fixed. Only
    the assignment of log ovule number to species is shuffled within families.
    The test asks whether the observed directional alignment is stronger than
    expected from family composition and the cross network alone.
    """
    data = pair_df.loc[pair_df["asymmetry_model_ready"]].copy().reset_index(drop=True)
    if data.empty or n_permutations <= 0:
        return pd.DataFrame({
            "test": ["family_stratified_ovule_permutation"],
            "n_permutations": [0],
            "observed_slope": [np.nan],
            "null_mean": [np.nan],
            "null_sd": [np.nan],
            "p_one_sided_positive": [np.nan],
        })

    q = build_control_q(data)
    y = data["sobel_chen_asymmetry_12"].to_numpy(dtype=float)
    y_res = residualize_against_controls(y, q)
    x_obs = data["delta_log_ovules_12"].to_numpy(dtype=float)
    x_obs_res = residualize_against_controls(x_obs, q)
    observed = float(np.dot(x_obs_res, y_res) / np.dot(x_obs_res, x_obs_res))

    rng = np.random.default_rng(seed)
    perm_slopes = np.empty(n_permutations, dtype=float)

    # Precompute family-specific species and original log ovules.
    family_maps: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for fam, g in data.groupby("family"):
        species_keys = pd.concat([g["species1_key"], g["species2_key"]]).dropna().unique()
        # First occurrence of a species has the same log ovule whichever side it appeared on.
        logs: dict[str, float] = {}
        for _, row in g.iterrows():
            logs[row["species1_key"]] = row["log_ovules_sp1"]
            logs[row["species2_key"]] = row["log_ovules_sp2"]
        species_keys = np.array([s for s in species_keys if s in logs], dtype=object)
        log_values = np.array([logs[s] for s in species_keys], dtype=float)
        family_maps[str(fam)] = (species_keys, log_values)

    for b in range(n_permutations):
        shuffled_logs: dict[str, float] = {}
        for fam, (species_keys, log_values) in family_maps.items():
            if len(species_keys) <= 1:
                for s, val in zip(species_keys, log_values):
                    shuffled_logs[s] = float(val)
            else:
                vals = rng.permutation(log_values)
                for s, val in zip(species_keys, vals):
                    shuffled_logs[s] = float(val)
        x_perm = np.array([
            shuffled_logs[row["species2_key"]] - shuffled_logs[row["species1_key"]]
            for _, row in data.iterrows()
        ], dtype=float)
        x_perm_res = residualize_against_controls(x_perm, q)
        denom = float(np.dot(x_perm_res, x_perm_res))
        perm_slopes[b] = np.nan if denom <= 0 else float(np.dot(x_perm_res, y_res) / denom)

    valid = perm_slopes[np.isfinite(perm_slopes)]
    p_one = (1.0 + np.sum(valid >= observed)) / (len(valid) + 1.0) if len(valid) else np.nan
    out = pd.DataFrame({
        "test": ["family_stratified_ovule_permutation"],
        "n_permutations": [int(len(valid))],
        "observed_slope_family_measure_residualized": [observed],
        "null_mean": [float(np.mean(valid)) if len(valid) else np.nan],
        "null_sd": [float(np.std(valid, ddof=1)) if len(valid) > 1 else np.nan],
        "p_one_sided_positive": [float(p_one) if np.isfinite(p_one) else np.nan],
        "null_slope_quantile_025": [float(np.quantile(valid, 0.025)) if len(valid) else np.nan],
        "null_slope_quantile_500": [float(np.quantile(valid, 0.500)) if len(valid) else np.nan],
        "null_slope_quantile_975": [float(np.quantile(valid, 0.975)) if len(valid) else np.nan],
    })
    return out


# -----------------------------------------------------------------------------
# Figures
# -----------------------------------------------------------------------------


def plot_funnel(funnel: pd.DataFrame, fig_dir: Path) -> None:
    plt.figure(figsize=(8, 5))
    y = np.arange(len(funnel))
    plt.barh(y, funnel["pair_rows"])
    plt.yticks(y, funnel["condition"])
    plt.xlabel("Pair rows")
    plt.gca().invert_yaxis()
    for yi, val in enumerate(funnel["pair_rows"]):
        plt.text(val + max(funnel["pair_rows"].max() * 0.01, 1), yi, str(int(val)), va="center", fontsize=8)
    plt.title("Data funnel for reciprocal asymmetry tests")
    save_current_figure(fig_dir, "data_funnel_reciprocal_asymmetry")


def plot_asymmetry_scatter(pair_df: pd.DataFrame, fig_dir: Path) -> None:
    d = pair_df.loc[pair_df["asymmetry_model_ready"]].copy()
    plt.figure(figsize=(7, 5))
    # Draw large families separately by default matplotlib colours; small families use the remaining default cycle.
    top_families = d["family"].value_counts().head(6).index.tolist()
    for fam in top_families:
        g = d.loc[d["family"] == fam]
        plt.scatter(g["delta_log_ovules_12"], g["sobel_chen_asymmetry_12"], s=28, alpha=0.75, label=fam)
    other = d.loc[~d["family"].isin(top_families)]
    if not other.empty:
        plt.scatter(other["delta_log_ovules_12"], other["sobel_chen_asymmetry_12"], s=20, alpha=0.35, label="Other families")
    plt.axhline(0, linewidth=0.8)
    plt.axvline(0, linewidth=0.8)
    plt.xlabel(r"$\Delta \log O_{12} = \log(O_2)-\log(O_1)$")
    plt.ylabel(r"$A^{SC}_{12}=RI_{1\leftarrow2}-RI_{2\leftarrow1}$")
    plt.title("Raw reciprocal asymmetry versus signed ovule contrast")
    plt.legend(fontsize=7, frameon=False, loc="best")
    save_current_figure(fig_dir, "asymmetry_vs_signed_ovule_contrast")


def plot_family_counts(family_summary: pd.DataFrame, fig_dir: Path, top_n: int = 20) -> None:
    f = family_summary.head(top_n).copy().iloc[::-1]
    plt.figure(figsize=(8, 6))
    y = np.arange(len(f))
    plt.barh(y, f["asymmetry_ready_pairs"])
    plt.yticks(y, f["family"])
    plt.xlabel("Asymmetry-ready pair rows")
    plt.title(f"Top {top_n} families for reciprocal asymmetry")
    save_current_figure(fig_dir, "family_asymmetry_ready_counts_top20")


def plot_family_slopes(family_slopes: pd.DataFrame, fig_dir: Path) -> None:
    if family_slopes.empty:
        return
    f = family_slopes.copy().sort_values("slope")
    f = f.loc[f["slope"].notna()].copy()
    plt.figure(figsize=(8, max(4, 0.28 * len(f))))
    y = np.arange(len(f))
    plt.errorbar(f["slope"], y, xerr=1.96 * f["standard_error"], fmt="o", markersize=4)
    plt.axvline(0, linewidth=0.8)
    plt.yticks(y, f["family"])
    plt.xlabel(r"Within-family slope of $A^{SC}$ on $\Delta\log O$")
    plt.title("Descriptive within-family slopes")
    save_current_figure(fig_dir, "descriptive_family_slopes")


def plot_leave_one_family_out(loo: pd.DataFrame, fig_dir: Path) -> None:
    if loo.empty:
        return
    d = loo.loc[loo["coefficient"].notna()].copy().sort_values("coefficient")
    plt.figure(figsize=(8, max(4, 0.22 * len(d))))
    y = np.arange(len(d))
    plt.errorbar(d["coefficient"], y, xerr=1.96 * d["standard_error"], fmt="o", markersize=3)
    plt.axvline(0, linewidth=0.8)
    plt.yticks(y, d["dropped_family"], fontsize=7)
    plt.xlabel(r"Slope after dropping one family")
    plt.title("Leave-one-family-out sensitivity")
    save_current_figure(fig_dir, "leave_one_family_out_slope")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyse reciprocal RI asymmetry for the ovule-number model.")
    parser.add_argument("--input-dir", type=Path, default=Path("/mnt/data"), help="Directory containing input CSV/tree files.")
    parser.add_argument("--output-dir", type=Path, default=Path("/mnt/data/ovule_prediction_analysis_update"), help="Directory for outputs.")
    parser.add_argument("--n-permutations", type=int, default=999, help="Number of within-family ovule-label permutations.")
    parser.add_argument("--use-exact-table2b-fill", action="store_true", help="Use exact non-removed species-level Table2b ovule values to fill missing Table2 ovules.")
    parser.add_argument("--seed", type=int, default=20260609, help="Random seed for permutation tests.")
    args = parser.parse_args()

    paths = Paths(args.input_dir, args.output_dir)
    paths.make_dirs()

    t1, t2, t2b, gd, tree = load_inputs(paths)
    all_table1_species = pd.concat([t1["species1"], t1["species2"]]).dropna().unique()
    ovule_lookup, fill_candidates = build_ovule_lookup(t2, t2b, all_table1_species)
    gd_collapsed, gd_audit = collapse_genetic_distances(gd)

    pair_df = make_pair_level_table(
        t1=t1,
        ovule_lookup=ovule_lookup,
        gd_collapsed=gd_collapsed,
        use_filled_ovules=args.use_exact_table2b_fill,
    )
    directional_df = make_directional_table(pair_df)

    funnel = make_funnel_table(pair_df)
    family_summary = make_family_summary(pair_df)
    family_modelability = make_family_modelability_table(family_summary)
    tree_audit = audit_tree(tree, pair_df)
    sign_tests = make_sign_tests(pair_df)
    model_summary = run_model_suite(pair_df, directional_df)
    family_slopes = make_family_slopes(pair_df)
    leave_one = make_leave_one_family_out(pair_df)
    permutation = family_stratified_ovule_permutation(pair_df, args.n_permutations, args.seed)

    # Negative-control summary: pairs with no ovule contrast.
    zero = pair_df.loc[pair_df["asymmetry_model_ready"] & (~pair_df["nonzero_ovule_difference"])].copy()
    nonzero = pair_df.loc[pair_df["asymmetry_model_ready"] & pair_df["nonzero_ovule_difference"]].copy()
    negative_control = pd.DataFrame({
        "subset": ["zero_ovule_contrast", "nonzero_ovule_contrast"],
        "n_pairs": [int(len(zero)), int(len(nonzero))],
        "mean_asymmetry": [float(zero["sobel_chen_asymmetry_12"].mean()) if len(zero) else np.nan, float(nonzero["sobel_chen_asymmetry_12"].mean()) if len(nonzero) else np.nan],
        "sd_asymmetry": [float(zero["sobel_chen_asymmetry_12"].std()) if len(zero) else np.nan, float(nonzero["sobel_chen_asymmetry_12"].std()) if len(nonzero) else np.nan],
        "mean_abs_asymmetry": [float(zero["sobel_chen_asymmetry_12"].abs().mean()) if len(zero) else np.nan, float(nonzero["sobel_chen_asymmetry_12"].abs().mean()) if len(nonzero) else np.nan],
    })

    # Write data.
    pair_df.to_csv(paths.data_dir / "reciprocal_pair_level_SC.csv", index=False)
    directional_df.to_csv(paths.data_dir / "directional_pair_level_SC.csv", index=False)
    ovule_lookup.to_csv(paths.data_dir / "lineage_ovule_lookup_with_provenance_fill_audit.csv", index=False)
    fill_candidates.to_csv(paths.table_dir / "exact_table2b_ovule_fill_candidates.csv", index=False)
    gd_collapsed.to_csv(paths.table_dir / "genetic_distance_collapsed_by_unordered_pair.csv", index=False)

    # Write summaries.
    funnel.to_csv(paths.table_dir / "data_funnel.csv", index=False)
    family_summary.to_csv(paths.table_dir / "family_summary.csv", index=False)
    family_modelability.to_csv(paths.table_dir / "family_modelability.csv", index=False)
    gd_audit.to_csv(paths.table_dir / "genetic_distance_audit.csv", index=False)
    tree_audit.to_csv(paths.table_dir / "tree_audit.csv", index=False)
    sign_tests.to_csv(paths.table_dir / "sign_and_correlation_tests.csv", index=False)
    model_summary.to_csv(paths.table_dir / "model_summary.csv", index=False)
    family_slopes.to_csv(paths.table_dir / "family_slopes.csv", index=False)
    leave_one.to_csv(paths.table_dir / "leave_one_family_out.csv", index=False)
    permutation.to_csv(paths.table_dir / "permutation_test.csv", index=False)
    negative_control.to_csv(paths.table_dir / "negative_control_zero_ovule_contrast.csv", index=False)

    # Figures.
    plot_funnel(funnel, paths.fig_dir)
    plot_asymmetry_scatter(pair_df, paths.fig_dir)
    plot_family_counts(family_summary, paths.fig_dir)
    plot_family_slopes(family_slopes, paths.fig_dir)
    plot_leave_one_family_out(leave_one, paths.fig_dir)

    # Plain-text note.
    ready_n = int(pair_df["asymmetry_model_ready"].sum())
    high_n = int((pair_df["asymmetry_model_ready"] & pair_df["both_confidence_le_2"]).sum())
    strict_n = int(pair_df["strict_PI_K2P_asymmetry_ready"].sum())
    readme = f"""Ovule-number reciprocal-asymmetry analysis outputs\n\nSign convention\n---------------\nA_12 = RI_SC(species 1 mother, species 2 father) - RI_SC(species 2 mother, species 1 father)\nDelta log O_12 = log(O_2) - log(O_1)\n\nA positive slope supports the prediction that the lower-ovule maternal species has the stronger prezygotic barrier.\n\nKey sample sizes\n----------------\nAsymmetry-ready reciprocal pairs: {ready_n}\nHigh-confidence ovule subset: {high_n}\nStrict PI/K2P subset: {strict_n}\n\nTable2b exact species-level fill used: {args.use_exact_table2b_fill}\nExact Table2b fill candidates found: {len(fill_candidates)}\n\nMain outputs\n------------\ndata/reciprocal_pair_level_SC.csv: one row per pair with exact Sobel-Chen directional RI and asymmetry.\ndata/directional_pair_level_SC.csv: two rows per pair for pair-fixed-effect models.\ntables/model_summary.csv: first-pass slope tests.\ntables/family_slopes.csv: descriptive within-family slopes.\ntables/leave_one_family_out.csv: leverage by family.\ntables/permutation_test.csv: within-family ovule-label permutation test.\nfigures/asymmetry_vs_signed_ovule_contrast.png: raw sign diagnostic.\n\nCaution\n-------\nThe first-pass models are checks, not the final paper model. The final model should use partial pooling, repeated-species effects, and phylogenetic covariance where possible.\n"""
    (paths.output_dir / "README.md").write_text(readme)

    print(readme)


if __name__ == "__main__":
    main()
