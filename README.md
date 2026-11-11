# Flat optics portfolio: FDTD, RCWA and GDS

Two small, fully reproducible projects that take thin-film photonics, the topic of my PhD, towards metasurface design.

| Project | What it shows | Tools |
|---|---|---|
| [01 – Microcavity: FDTD vs TMM](01_cavity_fdtd_vs_tmm) | A 1D FDTD model of an Ag/dielectric/Ag Fabry–Pérot cavity, validated against my own transfer-matrix code, with a grid-convergence study | Meep (FDTD), Python TMM |
| [02 – TiO₂ metalens: RCWA to GDS](02_metalens_rcwa_gds) | An RCWA unit-cell library for TiO₂ nanopillars, a 40 µm metalens designed from it, a focusing check, and a GDS layout ready for e-beam lithography | grcwa (RCWA), gdsfactory, KLayout |

Each project checks itself against a method it should agree with exactly:

- **FDTD vs TMM:** the cavity resonance agrees to within 0.01 nm. The RMS difference in reflectance drops from 4×10⁻³ to 3×10⁻⁴ as the grid is refined.
- **RCWA vs TMM:** for a uniform slab, R, T and the complex transmission coefficient agree to 6 decimal places.

<p align="center">
  <img src="01_cavity_fdtd_vs_tmm/results/fdtd_vs_tmm.png" width="45%">
  <img src="02_metalens_rcwa_gds/results/metalens_layout.png" width="52%">
</p>

## Quick start

```bash
# Meep is distributed through conda-forge
conda env create -f environment.yml
conda activate flat-optics
./run_all.sh            # about 3 minutes on a laptop
```

All figures and tables in `*/results/` are produced by these scripts.

## Repository layout

```
01_cavity_fdtd_vs_tmm/
    src/stack.py         cavity definition shared by both methods
    src/tmm.py           transfer-matrix method (s/p, oblique, absorbing media)
    src/run_fdtd.py      Meep 1D FDTD, two-run flux normalisation
    src/compare.py       overlay, Q-factor extraction, convergence, dispersion map
02_metalens_rcwa_gds/
    src/design.py              wavelength, materials, lattice and lens parameters
    src/unit_cell_library.py   RCWA sweep of pillar diameter + validation
    src/design_metalens.py     phase map, pillar choice, focusing check
    src/make_gds.py            GDS export with gdsfactory
```

## Author

Dr. Parvathi Valsalan, experimental physicist (PhD, Heidelberg University): thin-film microcavities, Fourier-space spectroscopy, ALD and thermal-evaporation fabrication.
