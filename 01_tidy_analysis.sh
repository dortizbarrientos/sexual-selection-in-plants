#!/usr/bin/env bash
set -uo pipefail

cd "$(git rev-parse --show-toplevel 2>/dev/null || echo .)"

if [ ! -d analysis ]; then
  echo "ERROR: run from the repo root (no analysis/ here)." >&2
  exit 1
fi

IS_GIT=0
git rev-parse --is-inside-work-tree >/dev/null 2>&1 && IS_GIT=1

gmv() {
  local src="$1" dst="$2"
  [ -e "$src" ] || { echo "  skip (missing): $src"; return 0; }
  if [ -e "$dst" ]; then echo "  skip (exists): $dst"; return 0; fi
  mkdir -p "$(dirname "$dst")"
  if [ $IS_GIT -eq 1 ] && [ -n "$(git ls-files "$src" 2>/dev/null)" ]; then
    git mv "$src" "$dst" && echo "  git mv: $src -> $dst"
  else
    mv "$src" "$dst" && echo "  mv    : $src -> $dst"
  fi
}

cd analysis

echo "==> [1] rename data audit: 01_2026_analyses -> 01_2026_data_audit"
gmv 01_2026_analyses 01_2026_data_audit

echo "==> [2] new self-contained reciprocal-asymmetry folder (script + outputs)"
mkdir -p 02_reciprocal_asymmetry_audit
gmv 03_ovule_project_audit 02_reciprocal_asymmetry_audit/outputs
gmv 04_reciprocal_family_asymmetry_audit/create_reciprocal_family_audit.py 02_reciprocal_asymmetry_audit/create_reciprocal_family_audit.py
rmdir 04_reciprocal_family_asymmetry_audit 2>/dev/null && echo "  removed empty: 04_reciprocal_family_asymmetry_audit"

echo "==> [3] home for the A_ij models: 02_2026_asymmetry -> 03_asymmetry_models"
rmdir 02_2026_asymmetry 2>/dev/null && echo "  removed empty: 02_2026_asymmetry"
mkdir -p 03_asymmetry_models
[ -e 03_asymmetry_models/.gitkeep ] || touch 03_asymmetry_models/.gitkeep

cd ..
echo "==> staging"
git add -A
git status
echo
echo "Next: git commit -m 'Tidy analysis/: data audit renamed, reciprocal-asymmetry script+outputs merged, asymmetry-models folder added'  &&  git push"
