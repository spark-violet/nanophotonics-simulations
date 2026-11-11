#!/usr/bin/env bash
# Reproduce every figure and table in the repository.
#   ./run_all.sh          projects 01-03 (03 takes ~1.5 h on 2 cores)
#   ./run_all.sh quick    projects 01-02 only (~3 min)
set -euo pipefail
cd "$(dirname "$0")"

echo "== 01: microcavity, FDTD vs TMM"
(cd 01_cavity_fdtd_vs_tmm/src && python run_fdtd.py && python compare.py)

echo "== 02: TiO2 metalens, RCWA -> design -> GDS"
(cd 02_metalens_rcwa_gds/src && python unit_cell_library.py && python design_metalens.py && python make_gds.py)

if [ "${1:-}" != "quick" ]; then
  echo "== 03: vortex-beam generator"
  (cd 03_vortex_beam_generator/src && python parameter_search.py && python pillar_library.py \
     && python fdtd_unit_cells.py && python inverse_design.py && python clean_freeform.py \
     && python vortex.py && python make_gds.py)
fi
echo "done"
