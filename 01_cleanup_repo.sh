#!/usr/bin/env bash
set -u

cd "$(git rev-parse --show-toplevel 2>/dev/null || echo .)"

echo "==> [1] Layer 1 engine -> src ; figures + grid -> results/demo/layer1"
if [ -d simulation/src/files ]; then
  [ -f simulation/src/files/layer1_competition_engine.py ] && mv -f simulation/src/files/layer1_competition_engine.py simulation/src/
  mkdir -p simulation/results/demo/layer1
  mv -f simulation/src/files/fig_layer1_*.png simulation/results/demo/layer1/ 2>/dev/null
  mv -f simulation/src/files/layer1_grid.csv   simulation/results/demo/layer1/ 2>/dev/null
  rmdir simulation/src/files 2>/dev/null
fi

echo "==> [2] remove download bundle + re-extracted duplicate package"
rm -f simulation/src/files.zip
DUP="simulation/src/ovule_gated_passenger_RI "
if [ -d "$DUP" ]; then
  if diff -rq "$DUP" archive/ovule_gated_passenger_RI >/dev/null 2>&1; then
    rm -rf "$DUP"; echo "    removed (identical to archive copy)"
  else
    echo "    KEPT: differs from archive/ovule_gated_passenger_RI -- inspect manually"
  fi
fi

echo "==> [3] retire superseded sim drafts -> archive/legacy_sims"
mkdir -p archive/legacy_sims
for f in sim1-june-8.py sim1-june-8.R sim2-june-8.py sim2-june-8.R simulate_ovule_gated_RI.py; do
  [ -f "simulation/src/$f" ] && mv -f "simulation/src/$f" archive/legacy_sims/
done

echo "==> [4] audit package -> analysis/2026_audit ; drop its zip"
if [ -d docs/ovule_project_audit ]; then
  mkdir -p analysis
  mv -f docs/ovule_project_audit analysis/2026_audit
fi
rm -f docs/ovule_project_audit_package.zip

echo "==> [5] verify pre-reorg twins are identical, then delete originals"
fail=0
chk() { if [ -e "$1" ] && [ -e "$2" ]; then diff -rq "$1" "$2" >/dev/null 2>&1 || { echo "    DIFFERS: $1  <>  $2"; fail=1; }; fi; }
for f in All_PI_pairs_Zuntini_2025_10.tre Genetic_distances_results.csv \
         Table1_2026_01_09.csv Table2_2025_12_15.csv Table2b_2025_11_25.csv \
         Table1_column_descriptions.txt Table2_column_descriptions.txt Table2b_column_descriptions.txt; do
  chk "2025_analyses/$f" "data/raw/2025_crossability/$f"
done
chk "2025_analyses/00_PI_ovule_number_analyses_2026_01_10.R" "analysis/2025_reproduction/00_PI_ovule_number_analyses_2026_01_10.R"
chk "2025_analyses/2025 Analyses Methods.docx" "analysis/2025_reproduction/2025 Analyses Methods.docx"
chk "2025_analyses/old" "archive/old"
chk "Ovule_RI_Task_List_June_2026.docx" "docs/Ovule_RI_Task_List_June_2026.docx"
chk "Reflective_Summary_Ovule_Gated_Model.tex" "theory/notes/Reflective_Summary_Ovule_Gated_Model.tex"
chk "sexual_selection_in_plants-notes.pdf" "theory/notes/sexual_selection_in_plants-notes.pdf"
chk "thoghts1.docx" "theory/notes/thoghts1.docx"

if [ "$fail" -eq 0 ]; then
  rm -rf 2025_analyses code
  rm -f Ovule_RI_Task_List_June_2026.docx Reflective_Summary_Ovule_Gated_Model.tex \
        sexual_selection_in_plants-notes.pdf thoghts1.docx
  echo "    all twins identical -- originals removed"
else
  echo "    SKIPPED deletion -- some twins differ (listed above)"
fi

echo "==> [6] tidy reorganize helper -> scripts/"
if [ -f 00_reorganize_repo.sh ]; then
  mkdir -p scripts
  mv -f 00_reorganize_repo.sh scripts/
fi

echo "==> staging"
git add -A
git status
echo
echo "Next: git commit -m 'Tidy repo: Layer 1 to src, outputs to results/demo, audit to analysis, remove duplicates'  &&  git push"
