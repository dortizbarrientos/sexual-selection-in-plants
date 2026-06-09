#!/usr/bin/env python3
"""Create reciprocal-cross family-level audit tables for the ovule-number project.

The script uses only the Python standard library for data handling so that it can
be run on any machine with Python 3.9+.
"""

from __future__ import annotations

import csv
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path('/mnt/data')
OUT = ROOT / 'ovule_project_audit'
OUT.mkdir(exist_ok=True)

TABLE1 = ROOT / 'Table1_2026_01_09.csv'
TABLE2 = ROOT / 'Table2_2025_12_15.csv'
GD = ROOT / 'Genetic_distances_results.csv'

MISSING = {'', 'NA', 'N/A', 'na', 'n/a', 'NaN', 'nan', 'None', 'none', 'NULL', 'null'}


def clean(x: object) -> str:
    return '' if x is None else str(x).strip()


def norm_name(x: object) -> str:
    return re.sub(r'\s+', ' ', clean(x).lower())


def parse_float(x: object) -> float | None:
    s = clean(x)
    if s in MISSING:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def pair_key(a: str, b: str) -> str:
    return ' || '.join(sorted([norm_name(a), norm_name(b)]))


def read_csv_dict(path: Path) -> list[dict[str, str]]:
    with path.open(newline='', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


# Load and index.
table1_all = read_csv_dict(TABLE1)
table1 = [r for r in table1_all if clean(r.get('remove')).lower() == 'no']

table2_rows = read_csv_dict(TABLE2)
ovules = {norm_name(r['species']): r for r in table2_rows}

gd_by_pair: dict[str, list[float | None]] = defaultdict(list)
for r in read_csv_dict(GD):
    gd_by_pair[pair_key(r['Species_1'], r['Species_2'])].append(parse_float(r['K2P_Distance']))


def is_pi(r: dict[str, str]) -> bool:
    return clean(r.get('PI_name')).upper() not in {'', 'NA', 'N/A'}


def has_reciprocal_inter(r: dict[str, str]) -> bool:
    return (
        parse_float(r.get('interspecific_value_sp1_sp2')) is not None
        and parse_float(r.get('interspecific_value_sp2_sp1')) is not None
    )


def has_directional_ri_inputs(r: dict[str, str]) -> bool:
    intra1 = parse_float(r.get('intraspecific_value_sp1'))
    intra2 = parse_float(r.get('intraspecific_value_sp2'))
    inter12 = parse_float(r.get('interspecific_value_sp1_sp2'))
    inter21 = parse_float(r.get('interspecific_value_sp2_sp1'))
    return (
        intra1 is not None and intra1 > 0
        and intra2 is not None and intra2 > 0
        and inter12 is not None
        and inter21 is not None
    )


def ovule_record(species: str) -> dict[str, str] | None:
    return ovules.get(norm_name(species))


def ovule_number(species: str) -> float | None:
    rec = ovule_record(species)
    if rec is None:
        return None
    v = parse_float(rec.get('ovule_number'))
    return v if v is not None and v > 0 else None


def confidence(species: str) -> float | None:
    rec = ovule_record(species)
    if rec is None:
        return None
    return parse_float(rec.get('confidence'))


def has_both_ovules(r: dict[str, str]) -> bool:
    return ovule_number(r['species1']) is not None and ovule_number(r['species2']) is not None


def has_high_confidence_ovules(r: dict[str, str]) -> bool:
    if not has_both_ovules(r):
        return False
    c1 = confidence(r['species1'])
    c2 = confidence(r['species2'])
    return c1 is not None and c2 is not None and c1 <= 2 and c2 <= 2


def k2p_values(r: dict[str, str]) -> list[float]:
    return [v for v in gd_by_pair.get(pair_key(r['species1'], r['species2']), []) if v is not None]


def has_k2p(r: dict[str, str]) -> bool:
    return len(k2p_values(r)) > 0


def directional_ri_values(r: dict[str, str]) -> tuple[float | None, float | None, float | None]:
    """Return RI species1-as-mother, RI species2-as-mother, and their difference.

    Uses column names as direction labels:
      interspecific_value_sp1_sp2 = species1 mother x species2 father
      interspecific_value_sp2_sp1 = species2 mother x species1 father
    """
    if not has_directional_ri_inputs(r):
        return None, None, None
    intra1 = parse_float(r['intraspecific_value_sp1'])
    intra2 = parse_float(r['intraspecific_value_sp2'])
    inter12 = parse_float(r['interspecific_value_sp1_sp2'])
    inter21 = parse_float(r['interspecific_value_sp2_sp1'])
    ri1 = 1.0 - inter12 / intra1
    ri2 = 1.0 - inter21 / intra2
    return ri1, ri2, ri1 - ri2


def signed_log_ovule_ratio(r: dict[str, str]) -> float | None:
    o1 = ovule_number(r['species1'])
    o2 = ovule_number(r['species2'])
    if o1 is None or o2 is None:
        return None
    return math.log10(o1) - math.log10(o2)


# Pair-level table, retaining only rows with reciprocal interspecific directions.
pair_rows: list[dict[str, object]] = []
for r in table1:
    if not has_reciprocal_inter(r):
        continue
    ri1, ri2, ri_asym = directional_ri_values(r)
    o1 = ovule_number(r['species1'])
    o2 = ovule_number(r['species2'])
    x = signed_log_ovule_ratio(r)
    gd_vals = k2p_values(r)
    pair_rows.append({
        'family': r['family'],
        'genus': r['genus'],
        'species1': r['species1'],
        'species2': r['species2'],
        'pair_key': pair_key(r['species1'], r['species2']),
        'PI_name': '' if clean(r.get('PI_name')).upper() == 'NA' else clean(r.get('PI_name')),
        'is_PI_from_PI_name': is_pi(r),
        'measure': r['measure'],
        'intraspecific_value_sp1': parse_float(r.get('intraspecific_value_sp1')),
        'intraspecific_value_sp2': parse_float(r.get('intraspecific_value_sp2')),
        'interspecific_value_sp1_sp2': parse_float(r.get('interspecific_value_sp1_sp2')),
        'interspecific_value_sp2_sp1': parse_float(r.get('interspecific_value_sp2_sp1')),
        'has_directional_RI_inputs': has_directional_ri_inputs(r),
        'directional_RI_sp1_mother': ri1,
        'directional_RI_sp2_mother': ri2,
        'directional_RI_asymmetry_sp1_minus_sp2': ri_asym,
        'ovule_sp1': o1,
        'ovule_sp2': o2,
        'confidence_sp1': confidence(r['species1']),
        'confidence_sp2': confidence(r['species2']),
        'has_both_ovules': has_both_ovules(r),
        'both_confidence_le_2': has_high_confidence_ovules(r),
        'signed_log10_ovule_ratio_sp1_minus_sp2': x,
        'abs_log10_ovule_difference': abs(x) if x is not None else None,
        'nonzero_ovule_difference': (abs(x) > 1e-12) if x is not None else False,
        'has_K2P': has_k2p(r),
        'K2P_Distance': gd_vals[0] if gd_vals else None,
        'asymmetry_model_ready': has_directional_ri_inputs(r) and has_both_ovules(r),
        'asymmetry_model_ready_confidence_le_2': has_directional_ri_inputs(r) and has_high_confidence_ovules(r),
        'strict_PI_K2P_asymmetry_ready': has_directional_ri_inputs(r) and has_both_ovules(r) and has_k2p(r),
    })

pair_fields = [
    'family', 'genus', 'species1', 'species2', 'pair_key', 'PI_name', 'is_PI_from_PI_name', 'measure',
    'intraspecific_value_sp1', 'intraspecific_value_sp2', 'interspecific_value_sp1_sp2', 'interspecific_value_sp2_sp1',
    'has_directional_RI_inputs', 'directional_RI_sp1_mother', 'directional_RI_sp2_mother', 'directional_RI_asymmetry_sp1_minus_sp2',
    'ovule_sp1', 'ovule_sp2', 'confidence_sp1', 'confidence_sp2', 'has_both_ovules', 'both_confidence_le_2',
    'signed_log10_ovule_ratio_sp1_minus_sp2', 'abs_log10_ovule_difference', 'nonzero_ovule_difference',
    'has_K2P', 'K2P_Distance', 'asymmetry_model_ready', 'asymmetry_model_ready_confidence_le_2', 'strict_PI_K2P_asymmetry_ready'
]
write_csv(OUT / 'reciprocal_pair_level_asymmetry_audit.csv', pair_rows, pair_fields)

# Family summary table.
family_names = sorted({r['family'] for r in table1})
family_rows: list[dict[str, object]] = []
for fam in family_names:
    fr = [r for r in table1 if r['family'] == fam]
    rec = [r for r in fr if has_reciprocal_inter(r)]
    directional = [r for r in rec if has_directional_ri_inputs(r)]
    asym = [r for r in directional if has_both_ovules(r)]
    high = [r for r in directional if has_high_confidence_ovules(r)]
    strict = [r for r in asym if has_k2p(r)]
    strict_high = [r for r in high if has_k2p(r)]

    species_recip = set()
    for r in rec:
        species_recip.update([r['species1'], r['species2']])
    species_asym = set()
    for r in asym:
        species_asym.update([r['species1'], r['species2']])

    x_vals = [signed_log_ovule_ratio(r) for r in asym if signed_log_ovule_ratio(r) is not None]
    x_high = [signed_log_ovule_ratio(r) for r in high if signed_log_ovule_ratio(r) is not None]
    nonzero = sum(1 for x in x_vals if abs(x) > 1e-12)
    nonzero_high = sum(1 for x in x_high if abs(x) > 1e-12)
    unique_x = len({round(x, 12) for x in x_vals})
    unique_x_high = len({round(x, 12) for x in x_high})

    genus_counts = Counter(r['genus'] for r in asym)
    top_genus, top_genus_n = ('', 0)
    if genus_counts:
        top_genus, top_genus_n = genus_counts.most_common(1)[0]

    n_species_asym = len(species_asym)
    possible_pairs = n_species_asym * (n_species_asym - 1) / 2 if n_species_asym >= 2 else 0
    density = len(asym) / possible_pairs if possible_pairs else None

    # Species reuse is a warning flag for pairwise non-independence.
    species_degree = Counter()
    for r in asym:
        species_degree[r['species1']] += 1
        species_degree[r['species2']] += 1
    max_species_reuse = max(species_degree.values()) if species_degree else 0

    family_rows.append({
        'family': fam,
        'retained_pairs': len(fr),
        'reciprocal_inter_pairs': len(rec),
        'directional_RI_pairs': len(directional),
        'asymmetry_pairs_with_ovules': len(asym),
        'asymmetry_pairs_with_nonzero_ovule_difference': nonzero,
        'unique_signed_log10_ovule_differences': unique_x,
        'high_confidence_asymmetry_pairs_conf_le_2': len(high),
        'high_confidence_nonzero_ovule_difference': nonzero_high,
        'unique_high_confidence_signed_log10_ovule_differences': unique_x_high,
        'strict_PI_K2P_asymmetry_pairs': len(strict),
        'strict_PI_K2P_high_confidence_pairs': len(strict_high),
        'species_in_reciprocal_pairs': len(species_recip),
        'species_in_asymmetry_ovule_pairs': n_species_asym,
        'genera_in_asymmetry_ovule_pairs': len(genus_counts),
        'top_genus_in_asymmetry_set': top_genus,
        'top_genus_asymmetry_pairs': top_genus_n,
        'top_genus_share_of_asymmetry_pairs': top_genus_n / len(asym) if asym else None,
        'asymmetry_pair_network_density_within_family_species': density,
        'max_times_one_species_reused_in_asymmetry_pairs': max_species_reuse,
        'family_has_at_least_2_asymmetry_pairs_and_2_ovule_differences': len(asym) >= 2 and unique_x >= 2,
        'family_has_at_least_3_asymmetry_pairs_and_2_ovule_differences': len(asym) >= 3 and unique_x >= 2,
        'family_has_at_least_5_asymmetry_pairs_and_2_ovule_differences': len(asym) >= 5 and unique_x >= 2,
    })

family_rows_sorted = sorted(
    family_rows,
    key=lambda r: (-int(r['asymmetry_pairs_with_ovules']), -int(r['reciprocal_inter_pairs']), str(r['family']))
)
family_fields = [
    'family', 'retained_pairs', 'reciprocal_inter_pairs', 'directional_RI_pairs',
    'asymmetry_pairs_with_ovules', 'asymmetry_pairs_with_nonzero_ovule_difference',
    'unique_signed_log10_ovule_differences', 'high_confidence_asymmetry_pairs_conf_le_2',
    'high_confidence_nonzero_ovule_difference', 'unique_high_confidence_signed_log10_ovule_differences',
    'strict_PI_K2P_asymmetry_pairs', 'strict_PI_K2P_high_confidence_pairs',
    'species_in_reciprocal_pairs', 'species_in_asymmetry_ovule_pairs', 'genera_in_asymmetry_ovule_pairs',
    'top_genus_in_asymmetry_set', 'top_genus_asymmetry_pairs', 'top_genus_share_of_asymmetry_pairs',
    'asymmetry_pair_network_density_within_family_species', 'max_times_one_species_reused_in_asymmetry_pairs',
    'family_has_at_least_2_asymmetry_pairs_and_2_ovule_differences',
    'family_has_at_least_3_asymmetry_pairs_and_2_ovule_differences',
    'family_has_at_least_5_asymmetry_pairs_and_2_ovule_differences',
]
write_csv(OUT / 'family_reciprocal_asymmetry_summary.csv', family_rows_sorted, family_fields)

# Global funnel table.
def count_family_threshold(rows: list[dict[str, object]], col: str, thresh: int) -> int:
    return sum(1 for r in rows if int(r[col]) >= thresh)

funnel = [
    {'stage': 'Retained pair rows', 'pairs': len(table1), 'families': len({r['family'] for r in table1})},
    {'stage': 'Both reciprocal interspecific directions', 'pairs': sum(has_reciprocal_inter(r) for r in table1), 'families': len({r['family'] for r in table1 if has_reciprocal_inter(r)})},
    {'stage': 'Directional RI possible: reciprocal + both intraspecific controls', 'pairs': sum(has_reciprocal_inter(r) and has_directional_ri_inputs(r) for r in table1), 'families': len({r['family'] for r in table1 if has_reciprocal_inter(r) and has_directional_ri_inputs(r)})},
    {'stage': 'Asymmetry model ready: directional RI + both ovule numbers', 'pairs': sum(has_directional_ri_inputs(r) and has_both_ovules(r) for r in table1), 'families': len({r['family'] for r in table1 if has_directional_ri_inputs(r) and has_both_ovules(r)})},
    {'stage': 'Asymmetry model ready with nonzero ovule difference', 'pairs': sum(has_directional_ri_inputs(r) and has_both_ovules(r) and abs(signed_log_ovule_ratio(r) or 0) > 1e-12 for r in table1), 'families': len({r['family'] for r in table1 if has_directional_ri_inputs(r) and has_both_ovules(r) and abs(signed_log_ovule_ratio(r) or 0) > 1e-12})},
    {'stage': 'High-confidence ovules only: confidence <= 2', 'pairs': sum(has_directional_ri_inputs(r) and has_high_confidence_ovules(r) for r in table1), 'families': len({r['family'] for r in table1 if has_directional_ri_inputs(r) and has_high_confidence_ovules(r)})},
    {'stage': 'Strict PI/K2P asymmetry subset', 'pairs': sum(has_directional_ri_inputs(r) and has_both_ovules(r) and has_k2p(r) for r in table1), 'families': len({r['family'] for r in table1 if has_directional_ri_inputs(r) and has_both_ovules(r) and has_k2p(r)})},
    {'stage': 'Strict PI/K2P + high-confidence ovules', 'pairs': sum(has_directional_ri_inputs(r) and has_high_confidence_ovules(r) and has_k2p(r) for r in table1), 'families': len({r['family'] for r in table1 if has_directional_ri_inputs(r) and has_high_confidence_ovules(r) and has_k2p(r)})},
]
write_csv(OUT / 'reciprocal_asymmetry_funnel.csv', funnel, ['stage', 'pairs', 'families'])

# Small human-readable text summary.
summary_path = OUT / 'reciprocal_asymmetry_readme.txt'
with summary_path.open('w', encoding='utf-8') as f:
    f.write('Ovule project reciprocal-cross family audit\n')
    f.write('============================================\n\n')
    f.write('Definitions used in this audit:\n')
    f.write('- Retained row: Table1 remove == no.\n')
    f.write('- Reciprocal interspecific pair: both interspecific_value_sp1_sp2 and interspecific_value_sp2_sp1 are numeric.\n')
    f.write('- Directional RI pair: reciprocal interspecific pair plus numeric positive intraspecific values for both species.\n')
    f.write('- Asymmetry model-ready pair: directional RI pair plus positive ovule_number for both species in Table2.\n')
    f.write('- High-confidence pair: asymmetry model-ready pair with both confidence values <= 2.\n')
    f.write('- Strict PI/K2P subset: asymmetry model-ready pair with nonmissing K2P distance. In these files this also maps to the PI subset.\n\n')
    for row in funnel:
        f.write(f"{row['stage']}: {row['pairs']} pairs across {row['families']} families\n")
    f.write('\nFamily-level modelability counts:\n')
    for label, col in [
        ('reciprocal interspecific pairs', 'reciprocal_inter_pairs'),
        ('directional RI pairs', 'directional_RI_pairs'),
        ('asymmetry pairs with ovules', 'asymmetry_pairs_with_ovules'),
        ('asymmetry pairs with nonzero ovule difference', 'asymmetry_pairs_with_nonzero_ovule_difference'),
        ('high-confidence asymmetry pairs', 'high_confidence_asymmetry_pairs_conf_le_2'),
        ('strict PI/K2P asymmetry pairs', 'strict_PI_K2P_asymmetry_pairs'),
    ]:
        f.write(f"\n{label}:\n")
        for thresh in [1, 2, 3, 5, 10, 20]:
            f.write(f"  families with >= {thresh}: {count_family_threshold(family_rows_sorted, col, thresh)}\n")

print('Wrote:')
for p in [
    OUT / 'family_reciprocal_asymmetry_summary.csv',
    OUT / 'reciprocal_pair_level_asymmetry_audit.csv',
    OUT / 'reciprocal_asymmetry_funnel.csv',
    summary_path,
]:
    print(p)
