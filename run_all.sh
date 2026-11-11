#!/usr/bin/env bash
# Reproduce every figure and table in the repository.
set -euo pipefail
cd "$(dirname "$0")"

echo "== 01: microcavity, FDTD vs TMM"
(cd 01_cavity_fdtd_vs_tmm/src && python run_fdtd.py && python compare.py)

echo "== 02: TiO2 metalens, RCWA -> design -> GDS"
(cd 02_metalens_rcwa_gds/src && python unit_cell_library.py && python design_metalens.py && python make_gds.py)

echo "done"
