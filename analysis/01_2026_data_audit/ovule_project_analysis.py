#!/usr/bin/env python3
"""
Ovule-number project: structural audit, model-readiness tables, and figures.

Inputs are expected in /mnt/data:
- Table1_2026_01_09.csv
- Table2_2025_12_15.csv
- Table2b_2025_11_25.csv
- Genetic_distances_results.csv
- All_PI_pairs_Zuntini_2025_10.tre
- crosses_template.csv
- lineage_traits_template.csv
- phylogeny_template.nwk

The script writes summary tables, figures, a LaTeX report, and model-ready joined
CSV files to /mnt/data/ovule_project_audit.
"""

from __future__ import annotations

import math
import os
import re
import textwrap
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

BASE = Path('/mnt/data')
OUT = BASE / 'ovule_project_audit'
FIG = OUT / 'figures'
TAB = OUT / 'tables'
REPORT = OUT / 'report'
for d in (OUT, FIG, TAB, REPORT):
    d.mkdir(parents=True, exist_ok=True)

# -----------------------------
# Helpers
# -----------------------------

def normalise_species_name(x: object) -> str | None:
    """A conservative species-name key used only for matching exact binomials."""
    if pd.isna(x):
        return None
    s = str(x).strip().replace('_', ' ')
    s = re.sub(r'\s+', ' ', s)
    return s.lower()


def pair_key(a: object, b: object) -> str:
    vals = [normalise_species_name(a), normalise_species_name(b)]
    return ' || '.join(sorted(vals))


def safe_divide(num: pd.Series, den: pd.Series) -> pd.Series:
    return num / den.replace({0: np.nan})


def savefig(name: str) -> None:
    for ext in ('pdf', 'png'):
        plt.savefig(FIG / f'{name}.{ext}', bbox_inches='tight', dpi=300)
    plt.close()


def percent(x: float) -> str:
    if pd.isna(x):
        return ''
    return f'{100*x:.1f}%'


def _format_value_for_latex(x: object) -> str:
    """Format values for compact LaTeX tables while leaving CSV files untouched."""
    if pd.isna(x):
        return ''
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    if isinstance(x, (float, np.floating)):
        if np.isfinite(x) and abs(x - round(x)) < 1e-9:
            return str(int(round(x)))
        return f'{float(x):.3f}'.rstrip('0').rstrip('.')
    return str(x)


def write_csv_and_tex(df: pd.DataFrame, stem: str, caption: str | None = None, label: str | None = None) -> None:
    df.to_csv(TAB / f'{stem}.csv', index=False)
    display_df = df.copy().applymap(_format_value_for_latex)
    latex = display_df.to_latex(index=False, escape=True, longtable=False, na_rep='', caption=caption, label=label)
    latex = latex.replace('\\begin{tabular}', '\\begin{adjustbox}{max width=\\textwidth}\n\\begin{tabular}', 1)
    latex = latex.replace('\\end{tabular}', '\\end{tabular}\n\\end{adjustbox}', 1)
    (TAB / f'{stem}.tex').write_text(latex)


def tex_escape(s: object) -> str:
    s = '' if pd.isna(s) else str(s)
    replacements = {
        '&': r'\&', '%': r'\%', '$': r'\$', '#': r'\#', '_': r'\_',
        '{': r'\{', '}': r'\}', '~': r'\textasciitilde{}', '^': r'\textasciicircum{}',
        '\\': r'\textbackslash{}'
    }
    out = ''
    for ch in s:
        out += replacements.get(ch, ch)
    return out

# -----------------------------
# Load data
# -----------------------------

t1 = pd.read_csv(BASE / 'Table1_2026_01_09.csv')
t2 = pd.read_csv(BASE / 'Table2_2025_12_15.csv')
t2b = pd.read_csv(BASE / 'Table2b_2025_11_25.csv')
gd = pd.read_csv(BASE / 'Genetic_distances_results.csv')
lineage_template = pd.read_csv(BASE / 'lineage_traits_template.csv')
crosses_template = pd.read_csv(BASE / 'crosses_template.csv')

# Clean flag fields without mutating original semantics.
t1['remove_flag'] = t1['remove'].fillna('no').str.lower().eq('yes')
t1['is_PI_from_PI_name'] = t1['PI_name'].notna()
t1['sp1_key'] = t1['species1'].map(normalise_species_name)
t1['sp2_key'] = t1['species2'].map(normalise_species_name)
t1['pair_key'] = [pair_key(a, b) for a, b in zip(t1['species1'], t1['species2'])]
t1_clean = t1.loc[~t1['remove_flag']].copy()
t1_pi = t1_clean.loc[t1_clean['is_PI_from_PI_name']].copy()

t2['species_key'] = t2['species'].map(normalise_species_name)
t2b['species_key'] = t2b['species'].map(normalise_species_name)
t2b['remove_flag'] = t2b['Remove'].fillna('no').str.lower().eq('yes')

gd['sp1_key'] = gd['Species_1'].map(normalise_species_name)
gd['sp2_key'] = gd['Species_2'].map(normalise_species_name)
gd['pair_key'] = [pair_key(a, b) for a, b in zip(gd['Species_1'], gd['Species_2'])]
gd_agg = gd.groupby('pair_key', as_index=False).agg(
    gd_rows=('K2P_Distance', 'size'),
    gd_nonmissing=('K2P_Distance', lambda x: int(x.notna().sum())),
    K2P_Distance=('K2P_Distance', 'mean'),
    K2P_min=('K2P_Distance', 'min'),
    K2P_max=('K2P_Distance', 'max')
)

# Lookup tables. Table2 species are unique.
ovule_lookup = t2.set_index('species_key')['ovule_number']
conf_lookup = t2.set_index('species_key')['confidence']
life_lookup = t2.set_index('species_key')['life_history']
mate_lookup = t2.set_index('species_key')['mating_system']
genus_lookup = t2.set_index('species_key')['genus']

# Enrich pair rows.
def enrich_pair_table(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out['ovule_sp1'] = out['sp1_key'].map(ovule_lookup)
    out['ovule_sp2'] = out['sp2_key'].map(ovule_lookup)
    out['confidence_sp1'] = out['sp1_key'].map(conf_lookup)
    out['confidence_sp2'] = out['sp2_key'].map(conf_lookup)
    out['life_history_sp1'] = out['sp1_key'].map(life_lookup)
    out['life_history_sp2'] = out['sp2_key'].map(life_lookup)
    out['mating_system_sp1'] = out['sp1_key'].map(mate_lookup)
    out['mating_system_sp2'] = out['sp2_key'].map(mate_lookup)
    out = out.merge(gd_agg, on='pair_key', how='left')
    out['has_both_ovule_numbers'] = out['ovule_sp1'].notna() & out['ovule_sp2'].notna()
    out['has_nonmissing_K2P'] = out['K2P_Distance'].notna()
    out['model_ready_PI_K2P_ovules'] = out['is_PI_from_PI_name'] & out['has_both_ovule_numbers'] & out['has_nonmissing_K2P'] & ~out['remove_flag']
    out['mean_log10_ovules'] = np.log10(np.sqrt(out['ovule_sp1'] * out['ovule_sp2']))
    out['abs_log10_ovule_difference'] = np.abs(np.log10(out['ovule_sp1']) - np.log10(out['ovule_sp2']))
    out['signed_log10_ovule_ratio_sp1_minus_sp2'] = np.log10(out['ovule_sp1']) - np.log10(out['ovule_sp2'])
    out['best_confidence_pair'] = out[['confidence_sp1', 'confidence_sp2']].max(axis=1)
    return out

pairs_all = enrich_pair_table(t1)
pairs_clean = pairs_all.loc[~pairs_all['remove_flag']].copy()
pairs_pi = pairs_clean.loc[pairs_clean['is_PI_from_PI_name']].copy()
model_ready_pairs = pairs_pi.loc[pairs_pi['model_ready_PI_K2P_ovules']].copy()
model_ready_pairs.to_csv(OUT / 'model_ready_PI_pairs.csv', index=False)
pairs_clean.to_csv(OUT / 'cleaned_pairs_with_ovules_and_distances.csv', index=False)

# Directional dataset skeleton: one row per maternal direction.
direction_rows = []
for _, r in pairs_clean.iterrows():
    for maternal_side in (1, 2):
        if maternal_side == 1:
            maternal, paternal = r['species1'], r['species2']
            maternal_ovules, paternal_ovules = r['ovule_sp1'], r['ovule_sp2']
            intra_maternal = r['intraspecific_value_sp1']
            inter_cross = r['interspecific_value_sp1_sp2']
            maternal_conf, paternal_conf = r['confidence_sp1'], r['confidence_sp2']
            maternal_life, paternal_life = r['life_history_sp1'], r['life_history_sp2']
            maternal_mating, paternal_mating = r['mating_system_sp1'], r['mating_system_sp2']
        else:
            maternal, paternal = r['species2'], r['species1']
            maternal_ovules, paternal_ovules = r['ovule_sp2'], r['ovule_sp1']
            intra_maternal = r['intraspecific_value_sp2']
            inter_cross = r['interspecific_value_sp2_sp1']
            maternal_conf, paternal_conf = r['confidence_sp2'], r['confidence_sp1']
            maternal_life, paternal_life = r['life_history_sp2'], r['life_history_sp1']
            maternal_mating, paternal_mating = r['mating_system_sp2'], r['mating_system_sp1']
        success_ratio = np.nan
        directional_RI = np.nan
        if pd.notna(inter_cross) and pd.notna(intra_maternal) and intra_maternal != 0:
            success_ratio = inter_cross / intra_maternal
            directional_RI = 1 - success_ratio
        direction_rows.append({
            'family': r['family'],
            'genus': r['genus'],
            'PI_name': r['PI_name'],
            'pair_key': r['pair_key'],
            'measure': r['measure'],
            'maternal_species': maternal,
            'paternal_species': paternal,
            'intra_maternal': intra_maternal,
            'inter_cross': inter_cross,
            'success_ratio_inter_over_intra': success_ratio,
            'directional_RI': directional_RI,
            'maternal_ovules': maternal_ovules,
            'paternal_ovules': paternal_ovules,
            'signed_log10_ovule_ratio_maternal_minus_paternal': (np.log10(maternal_ovules) - np.log10(paternal_ovules)) if pd.notna(maternal_ovules) and pd.notna(paternal_ovules) and maternal_ovules > 0 and paternal_ovules > 0 else np.nan,
            'maternal_confidence': maternal_conf,
            'paternal_confidence': paternal_conf,
            'maternal_life_history': maternal_life,
            'paternal_life_history': paternal_life,
            'maternal_mating_system': maternal_mating,
            'paternal_mating_system': paternal_mating,
            'K2P_Distance': r['K2P_Distance'],
            'is_PI_from_PI_name': r['is_PI_from_PI_name'],
        })

directional = pd.DataFrame(direction_rows)
directional['direction_complete_no_K2P'] = directional[['intra_maternal', 'inter_cross', 'maternal_ovules', 'paternal_ovules']].notna().all(axis=1)
directional['direction_complete_with_K2P'] = directional['direction_complete_no_K2P'] & directional['K2P_Distance'].notna()
directional.to_csv(OUT / 'directional_cleaned_crosses.csv', index=False)

# Parse tree tips.
tree_text = (BASE / 'All_PI_pairs_Zuntini_2025_10.tre').read_text().strip()
tips = re.findall(r'([A-Za-z][A-Za-z0-9_.-]*):[0-9.eE+-]+', tree_text)
tip_counts = Counter(tips)
tip_summary = pd.DataFrame({
    'metric': ['tip instances', 'unique tip labels', 'duplicated tip labels', 'PI names in clean Table1', 'PI names found in tree', 'PI names missing from tree'],
    'value': [len(tips), len(set(tips)), sum(1 for k, v in tip_counts.items() if v > 1), t1_pi['PI_name'].nunique(), len(set(t1_pi['PI_name'].dropna()) & set(tips)), len(set(t1_pi['PI_name'].dropna()) - set(tips))]
})
if any(v > 1 for v in tip_counts.values()):
    dup_tips = ', '.join(sorted([k for k, v in tip_counts.items() if v > 1]))
else:
    dup_tips = 'none'

# -----------------------------
# Summary tables
# -----------------------------

# Dataset inventory.
def unique_species_pair_df(df: pd.DataFrame) -> int:
    if {'species1', 'species2'} <= set(df.columns):
        return int(pd.concat([df['species1'], df['species2']]).nunique())
    if 'species' in df.columns:
        return int(df['species'].nunique())
    if {'Species_1', 'Species_2'} <= set(df.columns):
        return int(pd.concat([df['Species_1'], df['Species_2']]).nunique())
    return np.nan

inventory = pd.DataFrame([
    {'dataset': 'Table1_2026_01_09.csv', 'role': 'Pairwise cross outcomes and RI', 'rows': len(t1), 'columns': 15, 'families': t1['family'].nunique(), 'genera': t1['genus'].nunique(), 'unique_species': unique_species_pair_df(t1)},
    {'dataset': 'Table2_2025_12_15.csv', 'role': 'Curated ovule, life-history, and mating-system covariates', 'rows': len(t2), 'columns': 6, 'families': np.nan, 'genera': t2['genus'].nunique(), 'unique_species': unique_species_pair_df(t2)},
    {'dataset': 'Table2b_2025_11_25.csv', 'role': 'Ovule-number master table before final curation', 'rows': len(t2b), 'columns': 7, 'families': t2b['family'].nunique(), 'genera': t2b['genus'].nunique(), 'unique_species': unique_species_pair_df(t2b)},
    {'dataset': 'Genetic_distances_results.csv', 'role': 'K2P genetic distances for PI species pairs', 'rows': len(gd), 'columns': 3, 'families': np.nan, 'genera': np.nan, 'unique_species': unique_species_pair_df(gd)},
    {'dataset': 'All_PI_pairs_Zuntini_2025_10.tre', 'role': 'Phylogeny over PI-pair labels', 'rows': len(tips), 'columns': np.nan, 'families': np.nan, 'genera': np.nan, 'unique_species': len(set(tips))},
])
write_csv_and_tex(inventory, 'dataset_inventory', 'Input files and their modelling roles.', 'tab:dataset_inventory')

# Table 1 quality summary.
t1_quality = pd.DataFrame([
    {'quantity': 'All pair rows', 'count': len(t1), 'comment': 'Before applying remove flag'},
    {'quantity': 'Retained pair rows', 'count': len(t1_clean), 'comment': 'remove == no'},
    {'quantity': 'Rows flagged for removal', 'count': int(t1['remove_flag'].sum()), 'comment': 'Duplicates or lower-quality alternatives'},
    {'quantity': 'Rows with PI_name', 'count': int(t1['is_PI_from_PI_name'].sum()), 'comment': 'Available proxy for PI flag'},
    {'quantity': 'Retained PI rows', 'count': len(t1_pi), 'comment': 'remove == no and PI_name present'},
    {'quantity': 'Unique PI names retained', 'count': t1_pi['PI_name'].nunique(), 'comment': 'Should map to phylogeny tips'},
    {'quantity': 'Negative RI values', 'count': int((t1['RI'] < 0).sum()), 'comment': 'Interspecific performance can exceed intraspecific baseline'},
    {'quantity': 'RI values equal to 1', 'count': int((t1['RI'] == 1).sum()), 'comment': 'Boundary values affect beta regression'},
])
write_csv_and_tex(t1_quality, 'table1_quality_summary', 'Important structural features of Table 1.', 'tab:table1_quality')

measure_counts = pairs_clean.groupby('measure').agg(
    retained_rows=('RI', 'size'),
    PI_rows=('is_PI_from_PI_name', 'sum'),
    both_ovules=('has_both_ovule_numbers', 'sum'),
    nonmissing_K2P=('has_nonmissing_K2P', 'sum'),
    complete_PI_K2P_ovules=('model_ready_PI_K2P_ovules', 'sum'),
    RI_mean=('RI', 'mean'),
    RI_median=('RI', 'median'),
    RI_min=('RI', 'min'),
    RI_max=('RI', 'max'),
).reset_index()
for c in ['RI_mean', 'RI_median', 'RI_min', 'RI_max']:
    measure_counts[c] = measure_counts[c].round(3)
write_csv_and_tex(measure_counts, 'measure_readiness_summary', 'Readiness by cross-outcome measure.', 'tab:measure_readiness')

confidence_summary = t2.groupby('confidence').agg(
    species=('species', 'size'),
    ovule_numbers_present=('ovule_number', lambda x: int(x.notna().sum())),
    ovule_numbers_missing=('ovule_number', lambda x: int(x.isna().sum())),
    median_ovules=('ovule_number', 'median'),
    max_ovules=('ovule_number', 'max')
).reset_index().sort_values('confidence')
confidence_summary['median_ovules'] = confidence_summary['median_ovules'].round(2)
confidence_summary['max_ovules'] = confidence_summary['max_ovules'].round(2)
write_csv_and_tex(confidence_summary, 'ovule_confidence_summary', 'Ovule-number coverage by confidence category.', 'tab:confidence_summary')

life_mating = pd.crosstab(t2['life_history'], t2['mating_system']).reset_index()
write_csv_and_tex(life_mating, 'life_history_by_mating_system', 'Life history by mating system in Table 2.', 'tab:life_mating')

master_summary = pd.DataFrame([
    {'quantity': 'Rows in master table', 'count': len(t2b)},
    {'quantity': 'Rows marked table1 == yes', 'count': int((t2b['table1'] == 'yes').sum())},
    {'quantity': 'Rows marked for removal', 'count': int(t2b['remove_flag'].sum())},
    {'quantity': 'Nonremoved rows', 'count': int((~t2b['remove_flag']).sum())},
    {'quantity': 'Nonremoved rows with ovule number', 'count': int((~t2b['remove_flag'] & t2b['number'].notna()).sum())},
    {'quantity': 'Unique nonremoved species names', 'count': int(t2b.loc[~t2b['remove_flag'], 'species'].nunique())},
    {'quantity': 'Unique nonremoved table1 species names with number', 'count': int(t2b.loc[(~t2b['remove_flag']) & (t2b['table1'] == 'yes') & t2b['number'].notna(), 'species'].nunique())},
])
write_csv_and_tex(master_summary, 'table2b_master_summary', 'Summary of the ovule-number master table.', 'tab:master_summary')

join_rows = []
for label, df in [('All Table1 rows', pairs_all), ('Retained Table1 rows', pairs_clean), ('Retained PI rows', pairs_pi)]:
    join_rows.append({
        'subset': label,
        'rows': len(df),
        'both_species_in_Table2': int((df['sp1_key'].isin(set(t2['species_key'])) & df['sp2_key'].isin(set(t2['species_key']))).sum()),
        'both_ovules_numeric': int(df['has_both_ovule_numbers'].sum()),
        'K2P_row_present': int(df['gd_rows'].notna().sum()),
        'K2P_nonmissing': int(df['has_nonmissing_K2P'].sum()),
        'both_ovules_and_K2P': int((df['has_both_ovule_numbers'] & df['has_nonmissing_K2P']).sum()),
        'both_ovules_K2P_conf_le_3': int((df['has_both_ovule_numbers'] & df['has_nonmissing_K2P'] & df['confidence_sp1'].le(3) & df['confidence_sp2'].le(3)).sum()),
        'both_ovules_K2P_conf_le_2': int((df['has_both_ovule_numbers'] & df['has_nonmissing_K2P'] & df['confidence_sp1'].le(2) & df['confidence_sp2'].le(2)).sum()),
    })
join_coverage = pd.DataFrame(join_rows)
write_csv_and_tex(join_coverage, 'join_coverage_summary', 'Join coverage across cross, ovule, and genetic-distance data.', 'tab:join_coverage')

directional_summary = pd.DataFrame([
    {'subset': 'Retained Table1 directional rows', 'directions': len(directional),
     'complete_without_K2P': int(directional['direction_complete_no_K2P'].sum()),
     'complete_with_K2P': int(directional['direction_complete_with_K2P'].sum()),
     'complete_PI_with_K2P': int((directional['direction_complete_with_K2P'] & directional['is_PI_from_PI_name']).sum()),
     'negative_directional_RI': int((directional['directional_RI'] < 0).sum())},
])
write_csv_and_tex(directional_summary, 'directional_readiness_summary', 'Readiness of the reciprocal-direction data structure.', 'tab:directional_readiness')

# Missingness summary across main joined table.
missingness = []
for source_name, df in [('Table1', t1), ('Table2', t2), ('Table2b', t2b), ('Genetic distances', gd)]:
    for col in df.columns:
        if col.endswith('_key') or col in {'pair_key', 'remove_flag', 'is_PI_from_PI_name'}:
            continue
        missingness.append({
            'dataset': source_name,
            'column': col,
            'rows': len(df),
            'missing': int(df[col].isna().sum()),
            'missing_fraction': df[col].isna().mean()
        })
missingness = pd.DataFrame(missingness).sort_values(['dataset', 'missing_fraction'], ascending=[True, False])
missingness['missing_percent'] = (100 * missingness['missing_fraction']).round(1)
missingness_for_tex = missingness.loc[missingness['missing'] > 0, ['dataset', 'column', 'missing', 'missing_percent']].head(20)
write_csv_and_tex(missingness, 'missingness_all_columns', 'Missingness by dataset and column.', 'tab:missingness_all')
write_csv_and_tex(missingness_for_tex, 'missingness_top20', 'The 20 columns with nonzero missingness shown first by dataset.', 'tab:missingness_top20')

# Taxonomic coverage table.
tax_cov = pairs_clean.groupby('family').agg(
    retained_rows=('RI', 'size'),
    retained_PI_rows=('is_PI_from_PI_name', 'sum'),
    unique_genera=('genus', 'nunique'),
    unique_species=('species1', lambda x: np.nan) # placeholder, replaced below
).reset_index()
# Count unique species per family manually.
unique_species_by_family = []
for fam, sub in pairs_clean.groupby('family'):
    unique_species_by_family.append({'family': fam, 'unique_species': pd.concat([sub['species1'], sub['species2']]).nunique()})
unique_species_by_family = pd.DataFrame(unique_species_by_family)
tax_cov = tax_cov.drop(columns=['unique_species']).merge(unique_species_by_family, on='family')
tax_cov = tax_cov.sort_values('retained_rows', ascending=False)
write_csv_and_tex(tax_cov.head(20), 'top20_family_coverage', 'Top 20 families by retained pair rows.', 'tab:top20_families')
tax_cov.to_csv(TAB / 'family_coverage_all.csv', index=False)

# Ovule outliers and low-confidence high values.
ovule_outliers = t2.sort_values('ovule_number', ascending=False).head(15)[['species', 'ovule_number', 'confidence', 'life_history', 'mating_system']]
write_csv_and_tex(ovule_outliers, 'largest_ovule_numbers', 'Largest ovule-number estimates in Table 2.', 'tab:largest_ovules')

# Genetic-distance summary.
gd_summary = pd.DataFrame([
    {'quantity': 'Rows in genetic-distance table', 'count': len(gd), 'comment': ''},
    {'quantity': 'Unique unordered species pairs', 'count': gd['pair_key'].nunique(), 'comment': 'Duplicates are identical for eight pairs'},
    {'quantity': 'Rows with missing K2P', 'count': int(gd['K2P_Distance'].isna().sum()), 'comment': ''},
    {'quantity': 'Rows with nonmissing K2P', 'count': int(gd['K2P_Distance'].notna().sum()), 'comment': ''},
    {'quantity': 'Median nonmissing K2P', 'count': round(float(gd['K2P_Distance'].median()), 4), 'comment': ''},
    {'quantity': 'Maximum nonmissing K2P', 'count': round(float(gd['K2P_Distance'].max()), 4), 'comment': ''},
])
write_csv_and_tex(gd_summary, 'genetic_distance_summary', 'Summary of genetic-distance data.', 'tab:gd_summary')
write_csv_and_tex(tip_summary, 'phylogeny_tip_summary', 'Summary of PI-name matching to the phylogeny.', 'tab:tree_summary')

# Descriptive correlations in complete PI subset.
def corr_row(df: pd.DataFrame, x: str, y: str) -> dict:
    sub = df[[x, y]].replace([np.inf, -np.inf], np.nan).dropna()
    if len(sub) < 3:
        return {'x': x, 'y': y, 'n': len(sub), 'pearson': np.nan, 'spearman': np.nan}
    pearson = sub[x].corr(sub[y], method='pearson')
    spearman = sub[x].corr(sub[y], method='spearman')
    return {'x': x, 'y': y, 'n': len(sub), 'pearson': round(float(pearson), 3), 'spearman': round(float(spearman), 3)}
correlation_summary = pd.DataFrame([
    corr_row(model_ready_pairs, 'RI', 'K2P_Distance'),
    corr_row(model_ready_pairs, 'RI', 'mean_log10_ovules'),
    corr_row(model_ready_pairs, 'RI', 'abs_log10_ovule_difference'),
    corr_row(model_ready_pairs, 'K2P_Distance', 'abs_log10_ovule_difference'),
])
write_csv_and_tex(correlation_summary, 'descriptive_correlations_PI_complete', 'Descriptive correlations in the retained PI complete-case subset; these are diagnostics, not model tests.', 'tab:correlations')

# -----------------------------
# Figures
# -----------------------------

plt.rcParams.update({
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 10,
    'figure.titlesize': 13,
})

# Fig 1: Dataset map.
fig, ax = plt.subplots(figsize=(11, 6.6))
ax.set_axis_off()
boxes = {
    'Table 1\nCross outcomes + RI\n578 rows; 560 retained\n133 retained PI pairs': (0.05, 0.58, 0.26, 0.25),
    'Table 2\nCurated ovule traits\n408 species\n311 with ovule number': (0.37, 0.58, 0.25, 0.25),
    'Genetic distances\n141 rows; 133 unordered pairs\n110 retained PI pairs with K2P': (0.70, 0.58, 0.25, 0.25),
    'PI phylogeny\n134 tip instances\n133 unique PI labels\nCorymbia duplicated': (0.70, 0.20, 0.25, 0.23),
    'Model-ready PI subset\n87 pairs with RI, ovules, K2P\n57 with both confidence <= 2': (0.37, 0.20, 0.27, 0.23),
    'Prospective templates\nlineage traits\n+ reciprocal crosses\nneeded for count-level tests': (0.05, 0.20, 0.27, 0.23),
}
for text, (x, y, w, h) in boxes.items():
    rect = Rectangle((x, y), w, h, fill=False, linewidth=1.5)
    ax.add_patch(rect)
    ax.text(x + w / 2, y + h / 2, text, ha='center', va='center', wrap=True, fontsize=9)

def arrow(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle='->', mutation_scale=12, linewidth=1.2))
arrow(0.31, 0.70, 0.37, 0.70)
arrow(0.62, 0.70, 0.70, 0.70)
arrow(0.82, 0.58, 0.82, 0.43)
arrow(0.49, 0.58, 0.50, 0.43)
arrow(0.31, 0.32, 0.36, 0.32)
arrow(0.70, 0.32, 0.64, 0.32)
ax.set_title('Data structure and joins for the ovule-number project')
savefig('fig1_data_map')

# Fig 2: Join cascade.
cascade_labels = ['Retained\nTable 1', 'PI pairs', 'Both ovules', 'K2P present', 'Ovules + K2P', 'Conf. <= 2']
cascade_counts = [len(pairs_clean), len(pairs_pi), int(pairs_pi['has_both_ovule_numbers'].sum()), int(pairs_pi['has_nonmissing_K2P'].sum()), len(model_ready_pairs), int((model_ready_pairs['confidence_sp1'].le(2) & model_ready_pairs['confidence_sp2'].le(2)).sum())]
fig, ax = plt.subplots(figsize=(8.5, 5))
ax.bar(range(len(cascade_counts)), cascade_counts)
ax.set_xticks(range(len(cascade_counts)))
ax.set_xticklabels(cascade_labels)
ax.set_ylabel('Rows')
ax.set_title('Complete-case cascade for the PI model')
for i, v in enumerate(cascade_counts):
    ax.text(i, v + max(cascade_counts) * 0.02, str(v), ha='center', va='bottom')
ax.set_ylim(0, max(cascade_counts) * 1.18)
savefig('fig2_model_readiness_cascade')

# Fig 3: Missingness by key Table1 columns.
key_cols = ['intraspecific_value_sp1', 'intraspecific_value_sp2', 'interspecific_value_sp1_sp2', 'interspecific_value_sp2_sp1', 'sp1_hybrid_proportion', 'sp2_hybrid_proportion', 'RI', 'measure', 'PI_name']
miss_frac = t1[key_cols].isna().mean().sort_values()
fig, ax = plt.subplots(figsize=(8.5, 5.3))
ax.barh(range(len(miss_frac)), miss_frac.values * 100)
ax.set_yticks(range(len(miss_frac)))
ax.set_yticklabels([c.replace('_', ' ') for c in miss_frac.index])
ax.set_xlabel('Missing (%)')
ax.set_title('Missingness in key Table 1 columns')
for i, v in enumerate(miss_frac.values * 100):
    ax.text(v + 1, i, f'{v:.1f}%', va='center')
ax.set_xlim(0, 105)
savefig('fig3_table1_missingness')

# Fig 4: RI distributions by measure.
measures = list(pairs_clean['measure'].value_counts().index)
data = [pairs_clean.loc[pairs_clean['measure'] == m, 'RI'].dropna().values for m in measures]
fig, ax = plt.subplots(figsize=(9, 5.2))
ax.boxplot(data, labels=[textwrap.fill(m, 12) for m in measures], showmeans=True)
ax.axhline(0, linestyle='--', linewidth=1)
ax.set_ylabel('RI')
ax.set_title('RI distribution differs by measurement type')
ax.text(0.02, 0.02, 'Dashed line: RI = 0', transform=ax.transAxes, va='bottom')
savefig('fig4_RI_by_measure')

# Fig 5: Ovule distribution and confidence.
fig, ax = plt.subplots(figsize=(8.5, 5))
valid_ov = t2.loc[t2['ovule_number'].notna() & (t2['ovule_number'] > 0), 'ovule_number']
ax.hist(np.log10(valid_ov), bins=30)
ax.set_xlabel('log10(ovule number)')
ax.set_ylabel('Species')
ax.set_title('Ovule numbers span several orders of magnitude')
# Mark the largest two low-confidence Cattleya points.
for x in np.log10(t2.sort_values('ovule_number', ascending=False)['ovule_number'].dropna().head(2)):
    ax.axvline(x, linestyle=':', linewidth=1)
savefig('fig5_ovule_number_distribution')

# Fig 6: Table 2 confidence x data availability.
conf = confidence_summary.copy()
fig, ax = plt.subplots(figsize=(8, 5))
ax.bar(conf['confidence'].astype(str), conf['species'], label='species')
ax.bar(conf['confidence'].astype(str), conf['ovule_numbers_present'], label='with ovule number')
ax.set_xlabel('Confidence category')
ax.set_ylabel('Species')
ax.set_title('Ovule-number availability by confidence category')
ax.legend(frameon=False)
savefig('fig6_confidence_coverage')

# Fig 7: Genetic distance vs RI in complete PI subset.
fig, ax = plt.subplots(figsize=(7.5, 5.5))
for m, sub in model_ready_pairs.groupby('measure'):
    ax.scatter(sub['K2P_Distance'], sub['RI'], label=m, alpha=0.8)
ax.axhline(0, linestyle='--', linewidth=1)
ax.set_xlabel('K2P genetic distance')
ax.set_ylabel('RI')
ax.set_title('Complete PI pairs: genetic distance and RI')
ax.legend(frameon=False, fontsize=8)
savefig('fig7_K2P_vs_RI')

# Fig 8: RI vs ovule difference in complete PI subset.
fig, ax = plt.subplots(figsize=(7.5, 5.5))
for m, sub in model_ready_pairs.groupby('measure'):
    ax.scatter(sub['abs_log10_ovule_difference'], sub['RI'], label=m, alpha=0.8)
ax.axhline(0, linestyle='--', linewidth=1)
ax.set_xlabel('|log10 ovules species 1 - log10 ovules species 2|')
ax.set_ylabel('RI')
ax.set_title('Complete PI pairs: ovule-number contrast and RI')
ax.legend(frameon=False, fontsize=8)
savefig('fig8_ovule_difference_vs_RI')

# Fig 9: Family coverage top 15.
top_fam = tax_cov.head(15).sort_values('retained_rows')
fig, ax = plt.subplots(figsize=(8.5, 6.2))
ax.barh(range(len(top_fam)), top_fam['retained_rows'])
ax.set_yticks(range(len(top_fam)))
ax.set_yticklabels(top_fam['family'])
ax.set_xlabel('Retained Table 1 rows')
ax.set_title('Taxonomic concentration in retained cross data')
for i, v in enumerate(top_fam['retained_rows']):
    ax.text(v + 0.5, i, str(int(v)), va='center')
savefig('fig9_family_coverage')

# Fig 10: Directional RI vs maternal-paternal ovule log ratio.
dir_complete = directional.loc[directional['direction_complete_with_K2P'] & directional['directional_RI'].replace([np.inf, -np.inf], np.nan).notna()].copy()
fig, ax = plt.subplots(figsize=(7.5, 5.5))
for m, sub in dir_complete.groupby('measure'):
    ax.scatter(sub['signed_log10_ovule_ratio_maternal_minus_paternal'], sub['directional_RI'], label=m, alpha=0.8)
ax.axhline(0, linestyle='--', linewidth=1)
ax.axvline(0, linestyle=':', linewidth=1)
ax.set_xlabel('log10 maternal ovules - log10 paternal ovules')
ax.set_ylabel('Directional RI = 1 - interspecific / intraspecific')
ax.set_title('Directional structure for reciprocal cross models')
ax.legend(frameon=False, fontsize=8)
savefig('fig10_directional_RI_by_ovule_ratio')

# -----------------------------
# LaTeX report
# -----------------------------

complete_count = len(model_ready_pairs)
complete_conf2_count = int((model_ready_pairs['confidence_sp1'].le(2) & model_ready_pairs['confidence_sp2'].le(2)).sum())
neg_ri_count = int((t1['RI'] < 0).sum())
ri_min = t1['RI'].min()
ri_max = t1['RI'].max()
table2_ovule_present = int(t2['ovule_number'].notna().sum())
table2_ovule_missing = int(t2['ovule_number'].isna().sum())
gd_missing = int(gd['K2P_Distance'].isna().sum())

# Files relative to report directory; report will compile from OUT/report with fig/tables one level up.
report_tex = rf'''
\documentclass[11pt]{{article}}
\usepackage[a4paper,margin=2.2cm]{{geometry}}
\usepackage{{graphicx}}
\usepackage{{booktabs}}
\usepackage{{array}}
\usepackage{{longtable}}
\usepackage{{adjustbox}}
\usepackage{{float}}
\usepackage{{caption}}
\usepackage{{hyperref}}
\usepackage{{amsmath, amssymb}}
\usepackage{{enumitem}}
\usepackage{{xcolor}}
\setlength{{\parindent}}{{0pt}}
\setlength{{\parskip}}{{0.65em}}
\captionsetup{{font=small,labelfont=bf}}
\hypersetup{{colorlinks=true,linkcolor=black,urlcolor=blue,citecolor=black}}

\newcommand{{\conclusionbox}}[1]{{\vspace{{0.3em}}\noindent\fbox{{\begin{{minipage}}{{0.96\linewidth}}#1\end{{minipage}}}}\vspace{{0.4em}}}}

\title{{Ovule-number project: data structure, model-readiness, and analysis plan}}
\author{{Prepared for Daniel Ortiz-Barrientos}}
\date{{8 June 2026}}

\begin{{document}}
\maketitle

\section*{{Purpose}}

This report characterises the files supplied for the ovule-number project and asks a practical question: which parts of the data can test a model linking ovule number, reproductive isolation, genetic distance, life history, mating system, and phylogenetic structure? I treat Table 1 as the pairwise reproductive-isolation and cross-outcome table, Table 2 as the curated lineage-level covariate table, Table 2b as the master ovule-number database, the genetic-distance table as pair-level divergence information, and the Newick file as the phylogenetic structure for the phylogenetically independent pairs.

The working model assumed here is not a single fitted equation. It is the data structure needed to test whether ovule number, or contrasts in ovule number between species, predict cross success or reproductive isolation once genetic distance, mating system, life history, measurement type, and phylogenetic non-independence are controlled. This matters because the same biological idea can be tested in two different ways. A pair-level test asks whether species pairs with larger ovule-number contrasts have larger RI. A directional test asks whether the outcome differs when the many-ovule species is the maternal parent rather than the paternal parent. The second test is closer to the mechanism, but it needs reciprocal cross data on the right scale.

\conclusionbox{{\textbf{{Main result.}} The data are already close to a useful comparative analysis. After removing rows flagged for removal, there are {len(t1_clean)} retained Table 1 rows and {len(t1_pi)} retained rows with \texttt{{PI\_name}}. The most conservative complete-case subset with PI status, both ovule numbers, and nonmissing K2P distance contains {complete_count} pairs. If both ovule estimates are restricted to confidence categories 1--2, this decreases to {complete_conf2_count} pairs.}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.98\linewidth]{{../figures/fig1_data_map.pdf}}
\caption{{Overview of the supplied data files and the main joins needed for model testing.}}
\end{{figure}}

\section*{{Data inventory}}

\input{{../tables/dataset_inventory.tex}}

Table 1 contains the response variables: reciprocal intraspecific and interspecific cross outcomes, hybrid proportions for pollen-competition rows, RI, and measurement type. Table 2 contains the main lineage-level predictors: ovule number, confidence in the estimate, life history, and mating system. Table 2b is valuable as provenance and as a source for future curation, but Table 2 is the better analysis table because it supplies one best estimate per species.

A small naming issue should be fixed before formal analysis. The column-description file refers to a \texttt{{PI}} column, but the uploaded Table 1 does not contain a separate \texttt{{PI}} column. It contains \texttt{{PI\_name}}. I therefore treated a nonmissing \texttt{{PI\_name}} as the operational PI flag in this audit. That is reasonable for now because the retained PI rows and the tree labels agree, but I would add an explicit Boolean column called \texttt{{is\_PI}} before modelling.

\input{{../tables/table1_quality_summary.tex}}

\section*{{Coverage and model readiness}}

The most important structural result is the cascade from all retained cross rows to complete model-ready rows. Table 2 matches all retained Table 1 species by exact binomial names, but not all Table 2 rows have a numeric ovule estimate. The genetic-distance file covers the PI pairs, but {gd_missing} rows have missing K2P values. These two sources of missingness define the effective sample size for the strongest test.

\begin{{figure}}[H]
\centering
\includegraphics[width=0.82\linewidth]{{../figures/fig2_model_readiness_cascade.pdf}}
\caption{{The effective sample size after imposing phylogenetic independence, ovule-number availability, K2P availability, and high-confidence ovule estimates.}}
\end{{figure}}

\input{{../tables/join_coverage_summary.tex}}

The retained PI subset contains {len(pairs_pi)} rows. Of these, {int(pairs_pi['has_both_ovule_numbers'].sum())} have both ovule numbers and {int(pairs_pi['has_nonmissing_K2P'].sum())} have nonmissing K2P distances. The intersection is {complete_count} rows. This means a complete-case PI analysis is feasible, but it will have limited power for interactions among ovule number, genetic distance, mating system, life history, and measurement type. The model should therefore start simple and add complexity only when the residual structure justifies it.

\input{{../tables/measure_readiness_summary.tex}}

Measurement type is a major source of heterogeneity. Fruit set, seed set, seed number per fruit, and pollen competition are not the same response. They differ in the biological process measured and in the denominator. Fruit set measures whether a fruit forms. Seed set includes failed crosses when expressed as seeds per pollination. Seed number per fruit excludes failed crosses. Pollen competition uses hybrid proportions rather than direct seed or fruit counts. Pooling these four responses without a measurement-type term would mix mechanisms and can create regression slopes that are hard to interpret.

\begin{{figure}}[H]
\centering
\includegraphics[width=0.86\linewidth]{{../figures/fig4_RI_by_measure.pdf}}
\caption{{RI values by measurement type among retained Table 1 rows. The distributions differ enough that measurement type should be modelled, stratified, or treated as a sensitivity factor.}}
\end{{figure}}

\section*{{Ovule-number data}}

Table 2 contains {len(t2)} species. Numeric ovule estimates are present for {table2_ovule_present} species and missing for {table2_ovule_missing}. Confidence is useful because it gives a simple way to separate species-level estimates from weaker genus-level or indirect estimates. Most missing ovule numbers are in confidence categories 4 and 5.

\input{{../tables/ovule_confidence_summary.tex}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.82\linewidth]{{../figures/fig5_ovule_number_distribution.pdf}}
\caption{{Ovule numbers are highly skewed, so model terms based on ovule number should use a log scale. The two largest estimates are low-confidence values and should be checked in sensitivity analyses.}}
\end{{figure}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.78\linewidth]{{../figures/fig6_confidence_coverage.pdf}}
\caption{{Ovule-number availability by confidence category. Confidence categories 1--3 carry nearly all numeric ovule estimates.}}
\end{{figure}}

\input{{../tables/largest_ovule_numbers.tex}}

The main modelling consequence is simple. Raw ovule number should not enter the model linearly. It spans several orders of magnitude, and a few taxa can dominate raw-scale fits. Use \(\log_{{10}}(O)\), where \(O\) is ovule number, and then construct biologically interpretable contrasts such as
\[
\Delta O = \left| \log_{{10}}(O_1) - \log_{{10}}(O_2) \right|,
\]
for pair-level analyses, and
\[
D_{{m,p}} = \log_{{10}}(O_m) - \log_{{10}}(O_p),
\]
for directional analyses, where \(m\) and \(p\) denote maternal and paternal species.

\section*{{Genetic distance and phylogeny}}

\input{{../tables/genetic_distance_summary.tex}}

The genetic-distance table has {len(gd)} rows and {gd['pair_key'].nunique()} unique unordered species pairs. Eight species pairs are duplicated in the file; the duplicated distances are identical in this audit, so they do not change the values, but they should still be collapsed before modelling. The larger issue is missing distance values: {gd_missing} rows have no K2P value. Because the retained PI cross rows all have a corresponding distance-table row, the missingness is not a join problem. It is a missing-value problem.

\input{{../tables/phylogeny_tip_summary.tex}}

The tree is nearly aligned with the PI structure. All retained \texttt{{PI\_name}} labels are present in the tree. However, the tree contains {len(tips)} tip instances but only {len(set(tips))} unique tip labels; the duplicated tip label is \texttt{{{tex_escape(dup_tips)}}}. Many phylogenetic comparative tools require unique tip labels. Before fitting phylogenetic models, this duplicate label should be resolved, either by renaming the tips with unique suffixes or by checking whether one occurrence is unintended.

\section*{{What the current data can test}}

The pair-level test is ready for a first analysis. The clean dataset \texttt{{model\_ready\_PI\_pairs.csv}} contains {complete_count} rows with RI, both ovule numbers, and K2P distance. A useful first model would be
\[
RI_i = \alpha + \beta_1 \log_{{10}}(d_i + c) + \beta_2 \Delta O_i + \beta_3 \bar O_i + \gamma_{{measure(i)}} + \varepsilon_i,
\]
where \(d_i\) is genetic distance, \(c\) is a small constant if zero distances are retained, \(\Delta O_i\) is the absolute log ovule-number contrast, \(\bar O_i\) is the mean log ovule number of the pair, and \(\gamma_{{measure(i)}}\) is a measurement-type effect. I would not start with all life-history and mating-system interactions because the complete-case PI subset is too small for that many terms.

\begin{{figure}}[H]
\centering
\includegraphics[width=0.80\linewidth]{{../figures/fig7_K2P_vs_RI.pdf}}
\caption{{Exploratory relationship between K2P distance and RI in the complete PI subset. This is a diagnostic plot, not a fitted model.}}
\end{{figure}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.80\linewidth]{{../figures/fig8_ovule_difference_vs_RI.pdf}}
\caption{{Exploratory relationship between absolute log ovule-number contrast and RI in the complete PI subset.}}
\end{{figure}}

\input{{../tables/descriptive_correlations_PI_complete.tex}}

The directional test is biologically attractive because ovule number belongs to the maternal flower. For example, if the model predicts that a many-ovule maternal parent buffers cross failure, or alternatively exposes stronger post-pollination sorting, the reciprocal direction matters. The audit creates a file called \texttt{{directional\_cleaned\_crosses.csv}} with one row per maternal direction.

\input{{../tables/directional_readiness_summary.tex}}

The directional structure is promising but uneven. In the retained data there are {len(directional)} possible reciprocal directions, of which {int(directional['direction_complete_no_K2P'].sum())} have the interspecific outcome, maternal intraspecific baseline, and both ovule numbers. Only {int(directional['direction_complete_with_K2P'].sum())} also have K2P distance. Directional RI can be calculated as
\[
RI_{{m,p}} = 1 - \frac{{C_{{m,p}}}}{{C_{{m,m}}}},
\]
where \(C_{{m,p}}\) is the cross outcome with maternal species \(m\) and paternal species \(p\), and \(C_{{m,m}}\) is the intraspecific baseline for the maternal species. This definition is useful for diagnostics, but it inherits all problems of ratio data. It can be negative, it can exceed one if denominators are small or inconsistent, and it loses the original sampling variance.

\begin{{figure}}[H]
\centering
\includegraphics[width=0.80\linewidth]{{../figures/fig10_directional_RI_by_ovule_ratio.pdf}}
\caption{{Directional RI against the signed maternal-paternal ovule-number contrast. The vertical line marks equal ovule numbers between maternal and paternal species.}}
\end{{figure}}

\section*{{Strengths of the datasets}}

First, the project already separates response data from lineage-level covariates. That is good practice. Table 1 carries the pairwise RI and cross outcomes, while Table 2 carries ovule number, confidence, life history, and mating system. This separation makes joins explicit and allows uncertainty in ovule estimates to be modelled rather than hidden.

Second, the PI structure is coherent. The retained \texttt{{PI\_name}} values in Table 1 match the tree labels. This supports a conservative primary analysis based on the PI subset. It also makes it possible to compare three approaches: all retained rows, retained PI rows, and phylogenetic models using the full tree.

Third, reciprocal cross columns are present. This is important because ovule number is maternal in its immediate mechanism. A pair-level RI index can test whether species-pair contrasts matter, but reciprocal cross data can test whether the maternal side of the contrast matters.

Fourth, confidence categories are available for ovule number. This is useful for sensitivity analysis. A model fitted to all numeric ovule estimates can be compared with a stricter model using only confidence categories 1--2 or 1--3.

Fifth, the prospective templates are pointed in the right direction. The lineage template includes ovules, pollen production, pollen load, stigma size, style length, mating system, pollination syndrome, and flowering time. The cross template includes flower counts, ovules exposed, seeds matured, fruit set, pollen dose, pollen-tube information, F1 viability, F1 fertility, and genetic distance. Those fields would allow a much stronger count-level model than the current aggregated RI table.

\section*{{Weaknesses and consequences for the model}}

The main weakness is that the effective sample size drops when the biologically strongest conditions are imposed. A PI complete-case model with both ovule numbers and K2P distance has {complete_count} rows. This is sufficient for a simple analysis, but not for a saturated model with multiple interactions and random effects. A good first analysis should therefore test a small number of contrasts: genetic distance, absolute ovule-number difference, mean ovule number, and measurement type.

The second weakness is measurement heterogeneity. RI values derived from fruit set, seed set, seed number per fruit, and pollen competition do not share the same denominator. This affects the meaning of the slope. If ovule number predicts seed number per fruit but not fruit set, then pooling the responses may obscure the mechanism. Treat measurement type as a fixed effect, fit stratified models as sensitivity checks, or use a hierarchical model with measure-specific intercepts and slopes.

The third weakness is the loss of denominator information. Aggregated RI values hide the number of flowers, ovules, fruits, or seeds behind each estimate. This matters because an RI value based on 20 crosses should not carry the same precision as an RI value based on 500 crosses. The future cross-level template solves this. For the current data, use robust standard errors, nonparametric bootstrap over pairs, and sensitivity analyses by measurement type and confidence.

The fourth weakness is missing and uncertain ovule data. Table 2 has {table2_ovule_missing} missing ovule numbers, and some large estimates have low confidence. The practical consequence is that complete-case analyses may be biased if missingness is phylogenetically or biologically structured. Use confidence-weighted sensitivity analyses and consider multiple imputation, but do not let imputation create false precision. Report complete-case and imputed results side by side.

The fifth weakness is the response scale. The RI column ranges from {ri_min:.4f} to {ri_max:.4f}, with {neg_ri_count} negative values. Negative values are biologically meaningful because they indicate interspecific performance exceeding the intraspecific baseline. Boundary values at 1 are also present. A standard beta regression is therefore not valid on raw RI without transformation or censoring. Prefer models on log relative success where possible, or use Gaussian/robust models for RI with careful residual checks.

The sixth weakness is tree-label duplication. A duplicated \texttt{{Corymbia}} tip may break phylogenetic covariance construction. This is easy to fix, but it must be fixed before using tools such as \texttt{{ape}}, \texttt{{phylolm}}, \texttt{{MCMCglmm}}, \texttt{{brms}}, or \texttt{{RevBayes}}.

\section*{{Suggested analysis tools}}

For the current pair-level data, start with a transparent R or Python pipeline that creates a single analysis table. The minimum columns are species pair, \texttt{{PI\_name}}, RI, measure, both ovule numbers, ovule confidence categories, life history, mating system, and K2P distance. The file \texttt{{model\_ready\_PI\_pairs.csv}} created by this audit is a first version of that table.

For simple first-pass models, use linear models or robust linear models on RI, with \(\Delta O\), mean log ovule number, K2P distance, and measurement type as predictors. This gives a readable baseline. Follow this with permutation or bootstrap intervals over PI pairs. The bootstrap should resample pairs, not rows nested within the same pair.

For phylogenetic sensitivity, use phylogenetic generalised least squares. In R, \texttt{{ape}} can read and clean the Newick tree, \texttt{{caper}} or \texttt{{nlme}} can fit PGLS-style models, and \texttt{{phylolm}} can fit phylogenetic linear models. If you need measurement-error models or crossed random effects, use \texttt{{brms}} with a phylogenetic covariance matrix. The PI subset should remain the primary check because it is easy to explain.

For reciprocal cross models, reshape the data into maternal-direction rows. A useful diagnostic response is \(\log(C_{{m,p}} + c) - \log(C_{{m,m}} + c)\), where \(c\) is a small offset chosen before looking at the results. The key predictor is \(D_{{m,p}} = \log_{{10}}(O_m)-\log_{{10}}(O_p)\). Include pair as a clustering unit. If raw counts become available, replace the ratio response with a binomial, beta-binomial, Poisson, negative-binomial, or hurdle model depending on the denominator and outcome.

For missing ovule numbers, first report complete-case results. Then run sensitivity analyses using confidence thresholds: all numeric estimates, confidence 1--3, and confidence 1--2. Multiple imputation can be useful, but it should include genus, family, life history, mating system, and phylogenetic information. Imputation should not be used as the only analysis.

For future data collection, use the supplied templates but keep count denominators. The most useful fields are flowers crossed, ovules exposed, seeds matured, fruits formed, pollen dose, pollen tubes reaching the ovary, F1 viability, and F1 fertility. These allow the model to distinguish three stages: pollination success, seed maturation, and later hybrid fitness. That distinction is where the ovule-number model can become mechanistic rather than purely comparative.

\section*{{Recommended next steps}}

The immediate next step is to verify the PI flag. Add an explicit \texttt{{is\_PI}} column to Table 1 and keep \texttt{{PI\_name}} as the tree label. Then collapse duplicate genetic-distance rows and decide how to treat the {gd_missing} missing K2P values. Next, inspect the largest ovule-number estimates, especially low-confidence values, because they can dominate log-scale diagnostics less than raw-scale fits but still affect leverage.

After that, fit three models in order. First, fit a complete-case PI pair-level model with RI as the response. Second, repeat the model under stricter ovule-confidence thresholds. Third, fit a directional reciprocal-cross model using \texttt{{directional\_cleaned\_crosses.csv}}. The directional model is the most biologically informative, but the pair-level model is the safest first test.

\conclusionbox{{\textbf{{Bottom line.}} The datasets are strong enough for a careful first paper-level analysis, provided the claims are framed around a conservative PI complete-case test and sensitivity analyses. The strongest future version of the project would move from aggregated RI to count-level reciprocal crosses, because ovule number is a denominator and a maternal trait, not just a lineage-level covariate.}}

\appendix
\section*{{Additional structural tables}}

\input{{../tables/table2b_master_summary.tex}}

\input{{../tables/missingness_top20.tex}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.90\linewidth]{{../figures/fig3_table1_missingness.pdf}}
\caption{{Missingness in key Table 1 fields. Hybrid-proportion columns are sparse because they apply mainly to pollen-competition rows.}}
\end{{figure}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.86\linewidth]{{../figures/fig9_family_coverage.pdf}}
\caption{{Retained Table 1 rows are concentrated in a subset of plant families. This should be considered when interpreting broad comparative claims.}}
\end{{figure}}

\end{{document}}
'''

(REPORT / 'ovule_project_audit_report.tex').write_text(report_tex)

# Short README.
readme = f"""
Ovule-number project audit
==========================

Generated outputs
-----------------
- report/ovule_project_audit_report.tex: LaTeX source for the audit report.
- report/ovule_project_audit_report.pdf: compiled PDF, if pdflatex is available.
- figures/: PDF and PNG diagnostic figures.
- tables/: CSV and LaTeX versions of summary tables.
- model_ready_PI_pairs.csv: retained PI complete-case pair table with RI, ovules, and K2P.
- cleaned_pairs_with_ovules_and_distances.csv: retained pair table joined to Table2 and genetic distances.
- directional_cleaned_crosses.csv: one-row-per-maternal-direction table for reciprocal cross models.
- ovule_project_analysis.py: reproducible analysis script copied from /mnt/data.

Key audit result
----------------
Retained Table 1 rows: {len(t1_clean)}.
Retained PI rows from PI_name: {len(t1_pi)}.
Complete PI rows with both ovule numbers and nonmissing K2P: {complete_count}.
Complete PI rows with both ovule confidence categories <= 2: {complete_conf2_count}.
"""
(OUT / 'README.md').write_text(readme)

# Copy this script into output directory for self-contained archive.
this_script = Path(__file__)
if this_script.exists():
    (OUT / 'ovule_project_analysis.py').write_text(this_script.read_text())

print('Wrote outputs to', OUT)
print('Complete PI model-ready pairs:', complete_count)
print('Directional complete with K2P:', int(directional['direction_complete_with_K2P'].sum()))
