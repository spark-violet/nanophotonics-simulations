"""RCWA unit-cell library for a TiO2 nanopillar metalens (grcwa).

A metalens is built from sub-wavelength pillars. Each pillar diameter gives a
different transmitted phase. Light enters from the glass substrate and
leaves into air. This script computes, for one wavelength, the
complex transmission t = |t| exp(i*phi) of the zeroth diffraction order for a
square lattice of identical cylindrical pillars (local periodic approximation).

Checks included:
  1. A uniform TiO2 slab computed with RCWA is compared with the TMM from
     project 01 (R, T and the phase of t must agree).
  2. Convergence of the phase with the number of Fourier orders.

Outputs (in ../results):
  library.csv          diameter_nm, t_abs, phase_rad
  unit_cell_library.png
  rcwa_checks.md
"""
import sys
from pathlib import Path

import grcwa
import matplotlib.pyplot as plt
import numpy as np

import design

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2] / "01_cavity_fdtd_vs_tmm" / "src"))
import tmm  # noqa: E402

grcwa.set_backend("numpy")
RES = HERE.parents[1] / "results"
RES.mkdir(exist_ok=True)
NGRID = 128  # real-space grid of the pillar layer
# grcwa scales the incident plane wave differently from the usual Fresnel
# convention: t_Fresnel = n_incident * t_grcwa. Check 1 below verifies this
# against the TMM for a uniform slab.
AMP_SCALE = design.N_SUBSTRATE


def pillar_eps(diameter_um):
    """Permittivity grid of one unit cell with a centred circular pillar."""
    x = (np.arange(NGRID) + 0.5) / NGRID * design.PERIOD - design.PERIOD / 2
    X, Y = np.meshgrid(x, x, indexing="ij")
    inside = X**2 + Y**2 <= (diameter_um / 2) ** 2
    eps = np.where(inside, design.N_PILLAR**2, design.N_SURROUND**2)
    return eps.astype(float)


def solve(eps_grid=None, uniform_eps=None, n_orders=101):
    """Return (R, T, complex t of the zeroth order) for one unit cell."""
    P, wl = design.PERIOD, design.WAVELENGTH
    obj = grcwa.obj(n_orders, [P, 0], [0, P], 1 / wl, 0.0, 0.0, verbose=0)
    # Light comes from the glass substrate and leaves into air, as in a real
    # metalens that is illuminated through its substrate.
    obj.Add_LayerUniform(0.0, design.N_SUBSTRATE**2)    # incident medium (glass)
    if uniform_eps is not None:
        obj.Add_LayerUniform(design.HEIGHT, uniform_eps)
    else:
        obj.Add_LayerGrid(design.HEIGHT, NGRID, NGRID)
    obj.Add_LayerUniform(0.0, design.N_SURROUND**2)     # exit medium (air)
    obj.Init_Setup()
    if uniform_eps is None:
        obj.GridLayer_geteps(eps_grid.flatten())
    # x-polarised plane wave at normal incidence
    obj.MakeExcitationPlanewave(1, 0, 0, 0, order=0)
    R, T = obj.RT_Solve(normalize=1)
    E, _ = obj.Solve_FieldFourier(2, 0.0)
    i0 = next(k for k, g in enumerate(obj.G) if g[0] == 0 and g[1] == 0)
    return R, T, E[0][i0]


def main():
    lines = ["# RCWA checks", ""]

    # --- check 1: uniform slab vs TMM -----------------------------------
    R, T, t = solve(uniform_eps=design.N_PILLAR**2)
    n_list = [design.N_SUBSTRATE, design.N_PILLAR, design.N_SURROUND]
    R_t, T_t, _ = tmm.rt(n_list, [design.HEIGHT], design.WAVELENGTH)
    _, t_t = tmm.amplitudes(n_list, [design.HEIGHT], design.WAVELENGTH)
    lines += ["## 1. Uniform TiO2 slab: RCWA vs TMM", "",
              "| | R | T | abs(t) | phase(t) (rad) |", "|---|---|---|---|---|",
              f"| RCWA | {R:.6f} | {T:.6f} | {abs(AMP_SCALE*t):.6f} | {np.angle(t):.6f} |",
              f"| TMM  | {R_t:.6f} | {T_t:.6f} | {abs(t_t):.6f} | {np.angle(t_t):.6f} |", ""]

    # --- check 2: convergence in Fourier orders ---------------------------
    d_test = 0.5 * (design.D_MIN + design.D_MAX)
    eps = pillar_eps(d_test)
    lines += [f"## 2. Convergence with Fourier orders (pillar diameter {d_test*1e3:.0f} nm)",
              "", "| Fourier orders | T (all orders) | phase(t) (rad) |", "|---|---|---|"]
    for n in [21, 51, 101, 151, 201, 301]:
        _, T, t = solve(eps_grid=eps, n_orders=n)
        lines.append(f"| {n} | {T:.4f} | {np.angle(t):.4f} |")
    lines.append("")

    # --- library sweep ---------------------------------------------------
    diameters = np.round(np.arange(design.D_MIN, design.D_MAX + 1e-9, 0.002), 4)
    t_list = []
    for d in diameters:
        _, _, t = solve(eps_grid=pillar_eps(d), n_orders=design.N_ORDERS)
        t_list.append(t)
    t_arr = np.array(t_list)
    phase = np.unwrap(np.angle(t_arr))
    phase -= phase[0]  # relative to the thinnest pillar
    np.savetxt(RES / "library.csv",
               np.column_stack([diameters * 1e3, np.abs(t_arr), phase,
                                t_arr.real, t_arr.imag]),
               delimiter=",", header="diameter_nm,t_abs,phase_rad,t_re,t_im",
               comments="")
    # power transmittance of the zeroth order (see AMP_SCALE)
    T0 = design.N_SURROUND / design.N_SUBSTRATE * np.abs(AMP_SCALE * t_arr) ** 2
    span = (phase.max() - phase.min()) / (2 * np.pi)
    lines += ["## 3. Library", "",
              f"- {len(diameters)} diameters from {design.D_MIN*1e3:.0f} to "
              f"{design.D_MAX*1e3:.0f} nm",
              f"- phase coverage: {span:.2f} x 2π",
              f"- zeroth-order transmittance T0: mean "
              f"{np.mean(T0):.3f}, min {np.min(T0):.3f}", ""]
    (RES / "rcwa_checks.md").write_text("\n".join(lines))
    print("\n".join(lines))

    fig, ax = plt.subplots(figsize=(5.6, 3.6), dpi=150)
    ax.plot(diameters * 1e3, phase / (2 * np.pi), color="#7a1530", lw=2)
    ax.set_xlabel("Pillar diameter (nm)")
    ax.set_ylabel("Phase / 2π", color="#7a1530")
    ax2 = ax.twinx()
    ax2.plot(diameters * 1e3, T0, color="#1f4e79", lw=1.5)
    ax2.set_ylabel("Transmittance (zeroth order)", color="#1f4e79")
    ax2.set_ylim(0, 1.05)
    ax.set_title(f"TiO$_2$ pillars, H = {design.HEIGHT*1e3:.0f} nm, "
                 f"P = {design.PERIOD*1e3:.0f} nm, λ = {design.WAVELENGTH*1e3:.0f} nm",
                 fontsize=9)
    fig.tight_layout()
    fig.savefig(RES / "unit_cell_library.png")


if __name__ == "__main__":
    main()
