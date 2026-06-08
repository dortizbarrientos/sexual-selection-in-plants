#!/usr/bin/env bash
# =============================================================================
# 00_reorganize_repo.sh
# Restructure the `sexual-selection-plants` repository into a clean research
# compendium: data / theory / simulation / analysis / results / docs / archive.
#
# SAFETY GUARANTEES
#   * Never deletes anything (no `rm`). Superseded material is *moved* to archive/.
#   * Uses `git mv` when a file is tracked (preserves history); plain `mv` otherwise.
#   * Idempotent: re-running skips files already moved or already present.
#   * Cherry-picked reusable assets are *copied* (cp -n), leaving originals intact
#     until the source package is archived at the end.
#
# USAGE
#   cd /path/to/sexual-selection-plants
#   bash 00_reorganize_repo.sh
# =============================================================================
set -uo pipefail

# --- guard: must be run from the repo root -----------------------------------
if [[ ! -d "2025_analyses" && ! -d "code" && ! -d "ovule_gated_passenger_RI" ]]; then
  echo "ERROR: run this from the sexual-selection-plants/ repo root" >&2
  echo "       (none of 2025_analyses/, code/, ovule_gated_passenger_RI/ found here)." >&2
  exit 1
fi

IS_GIT=0
git rev-parse --is-inside-work-tree >/dev/null 2>&1 && IS_GIT=1

# --- helpers -----------------------------------------------------------------
# mv_safe <src> <dst-dir> : move src INTO dst-dir/, safely & idempotently
mv_safe() {
  local src="$1" dstdir="$2" base
  if [[ ! -e "$src" ]]; then echo "  skip (missing): $src"; return 0; fi
  base="$(basename "$src")"
  mkdir -p "$dstdir"
  if [[ -e "$dstdir/$base" ]]; then echo "  skip (exists) : $dstdir/$base"; return 0; fi
  if [[ $IS_GIT -eq 1 && -n "$(git ls-files "$src" 2>/dev/null)" ]]; then
    git mv "$src" "$dstdir/" && echo "  git mv: $src -> $dstdir/"
  else
    mv "$src" "$dstdir/"     && echo "  mv    : $src -> $dstdir/"
  fi
}

# cp_safe <src> <dst-dir> : copy src INTO dst-dir/ without overwriting (cp -n)
cp_safe() {
  local src="$1" dstdir="$2"
  if [[ ! -e "$src" ]]; then echo "  skip (missing): $src"; return 0; fi
  mkdir -p "$dstdir"
  cp -n "$src" "$dstdir/" && echo "  cp    : $src -> $dstdir/"
}

keep() { mkdir -p "$1"; [[ -e "$1/.gitkeep" ]] || touch "$1/.gitkeep"; }

echo "=== Reorganizing sexual-selection-plants (git mode: $IS_GIT) ==="

echo "[1/9] raw data  ->  data/raw/2025_crossability/"
for f in All_PI_pairs_Zuntini_2025_10.tre Genetic_distances_results.csv \
         Table1_2026_01_09.csv Table2_2025_12_15.csv Table2b_2025_11_25.csv \
         Table1_column_descriptions.txt Table2_column_descriptions.txt \
         Table2b_column_descriptions.txt; do
  mv_safe "2025_analyses/$f" "data/raw/2025_crossability"
done

echo "[2/9] 2025 reproduction analysis  ->  analysis/2025_reproduction/"
mv_safe "2025_analyses/00_PI_ovule_number_analyses_2026_01_10.R" "analysis/2025_reproduction"
mv_safe "2025_analyses/2025 Analyses Methods.docx"               "analysis/2025_reproduction"

echo "[3/9] theory documents  ->  theory/"
mv_safe "ovule_RI_theory_note.pdf"               "theory/theory_note"
mv_safe "Reflective_Summary_Ovule_Gated_Model.tex" "theory/notes"
mv_safe "thoghts1.docx"                          "theory/notes"
mv_safe "sexual_selection_in_plants-notes.pdf"   "theory/notes"
# manuscript skeleton + bibliography: COPIED from the autogen package (archived below)
cp_safe "ovule_gated_passenger_RI/latex/main.tex"        "theory/manuscript"
cp_safe "ovule_gated_passenger_RI/latex/references.bib"  "theory/manuscript"
cp_safe "ovule_gated_passenger_RI/latex/main.pdf"        "theory/manuscript"

echo "[4/9] simulation code  ->  simulation/src/"
for f in sim1-june-8.R sim1-june-8.py sim2-june-8.R sim2-june-8.py; do
  mv_safe "code/$f" "simulation/src"
done
# the autogen simulation engine: COPIED in for reference / consolidation
cp_safe "ovule_gated_passenger_RI/simulation/simulate_ovule_gated_RI.py" "simulation/src"

echo "[5/9] analysis R templates  ->  analysis/R/"
cp_safe "ovule_gated_passenger_RI/simulation/analysis_brms_template.R" "analysis/R"

echo "[6/9] data templates  ->  data/templates/"
for f in crosses_template.csv lineage_traits_template.csv phylogeny_template.nwk; do
  cp_safe "ovule_gated_passenger_RI/data_templates/$f" "data/templates"
done

echo "[7/9] canonical environment  ->  ./requirements.txt"
mv_safe "code/requirements.txt" "."   # NOTE: reconcile with archived package's requirements.txt

echo "[8/9] docs  ->  docs/"
mv_safe "Ovule_RI_Task_List_June_2026.docx" "docs"

echo "[9/9] archive superseded / unvetted material  ->  archive/"
mv_safe "2025_analyses/old"                  "archive"   # -> archive/old/
mv_safe "ovule_gated_passenger_RI"           "archive"   # whole autogen package, intact
mv_safe "ovule_gated_passenger_RI_package.zip" "archive"

echo "[+]   placeholders for empty dirs we will populate"
keep data/processed
keep simulation/configs
keep simulation/results
keep analysis/asymmetry
keep results/figures

echo
echo "=== Done. ==="
echo "Manual follow-ups (intentionally NOT automated):"
echo "  * Empty leftover dirs may remain (2025_analyses/, code/). Remove when happy:"
echo "      rmdir 2025_analyses code 2>/dev/null"
echo "  * Two requirements.txt existed (code/ and the autogen package). The first is now"
echo "      ./requirements.txt; the second is in archive/ovule_gated_passenger_RI/. Reconcile."
echo "  * Add the theory-note SOURCE (ovule_RI_theory_note.tex) to theory/theory_note/."
echo "  * Decide whether archive/ovule_gated_passenger_RI/ is the simulation base or a"
echo "      reference only. Its results_demo/ was generated wholesale and is unvetted."
