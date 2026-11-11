"""3D FDTD (Meep) check of the RCWA unit cells.

RCWA and FDTD solve Maxwell's equations in completely different ways
(Fourier modes in the plane vs. a time-stepped grid in real space). If both
give the same transmitted phase and amplitude for the 8 chosen pillars, the
library can be trusted.

Each unit cell is simulated in 3D with periodic boundaries in x and y and
PML in z. A broadband x-polarised plane wave is launched in the glass; the
transmitted field is recorded at the design wavelength on a plane in air
and averaged over the cell (= zeroth diffraction order).

  phase and |t|: relative to an empty cell (no pillar), in both methods
  T0          : transmitted flux / incident flux (homogeneous-glass run)

Usage:  python fdtd_unit_cells.py [resolution ...]   (default 40 80)
Outputs: ../results/fdtd_cells_res<R>.csv and fdtd_vs_rcwa.png
"""
import sys
from pathlib import Path

import grcwa
import matplotlib.pyplot as plt
import meep as mp
import numpy as np

import config as C
import rcwa
from pillar_library import load_geometry

mp.verbosity(0)
grcwa.set_backend("numpy")
RES = Path(__file__).resolve().parents[1] / "results"

PML, GLASS, AIR = 0.5, 0.5, 0.6
FCEN = 1 / C.WAVELENGTH


def run_cell(diameter, height, period, resolution, all_glass=False):
    """Return (complex Ex averaged over the monitor plane, transmitted flux)."""
    sz = 2 * PML + GLASS + height + AIR
    z_sub_top = -sz / 2 + PML + GLASS            # glass / pillar interface
    glass = mp.Medium(index=C.N_SUBSTRATE)
    geom = []
    if all_glass:
        default = glass
    else:
        default = mp.air
        geom.append(mp.Block(size=mp.Vector3(mp.inf, mp.inf, PML + GLASS),
                             center=mp.Vector3(z=-sz / 2 + (PML + GLASS) / 2),
                             material=glass))
        if diameter > 0:
            geom.append(mp.Cylinder(radius=diameter / 2, height=height,
                                    axis=mp.Vector3(z=1),
                                    center=mp.Vector3(z=z_sub_top + height / 2),
                                    material=mp.Medium(index=C.N_PILLAR)))
    z_src = -sz / 2 + PML + 0.15
    z_mon = z_sub_top + height + 0.35
    src = [mp.Source(mp.GaussianSource(FCEN, fwidth=0.3 * FCEN), component=mp.Ex,
                     center=mp.Vector3(z=z_src), size=mp.Vector3(period, period, 0))]
    sim = mp.Simulation(cell_size=mp.Vector3(period, period, sz),
                        boundary_layers=[mp.PML(PML, direction=mp.Z)],
                        geometry=geom, sources=src, default_material=default,
                        k_point=mp.Vector3(), resolution=resolution)
    plane = mp.Volume(center=mp.Vector3(z=z_mon), size=mp.Vector3(period, period, 0))
    dft = sim.add_dft_fields([mp.Ex], FCEN, 0, 1, where=plane)
    flux = sim.add_flux(FCEN, 0, 1, mp.FluxRegion(center=mp.Vector3(z=z_mon),
                                                   size=mp.Vector3(period, period, 0)))
    sim.run(until_after_sources=mp.stop_when_dft_decayed(1e-7))
    ex = np.mean(sim.get_dft_array(dft, mp.Ex, 0))
    return ex, mp.get_fluxes(flux)[0]


def main():
    H, P = load_geometry()
    levels = np.loadtxt(RES / "levels_pillars.csv", delimiter=",", skiprows=1)
    diams = levels[:, 2] / 1e3
    resolutions = [int(r) for r in sys.argv[1:]] or [40, 80]

    # RCWA reference, relative to the empty cell
    t_empty = rcwa.solve(rcwa.disk_eps(0.0, P, 128), H, P, 151)[1]
    t_rcwa = levels[:, 5] + 1j * levels[:, 6]
    rel_rcwa = t_rcwa / t_empty

    for res in resolutions:
        out = RES / f"fdtd_cells_res{res}.csv"
        if out.exists():
            print(f"res {res}: using existing {out.name} (delete it to rerun)")
            continue
        ex_empty, _ = run_cell(0.0, H, P, res)
        _, flux_inc = run_cell(0.0, H, P, res, all_glass=True)
        rows = []
        for k, d in enumerate(diams):
            ex, flux = run_cell(d, H, P, res)
            rel = ex / ex_empty
            rows.append((k, d * 1e3, flux / flux_inc, abs(rel), np.angle(rel),
                         rcwa.t_to_T0(t_rcwa[k]), abs(rel_rcwa[k]), np.angle(rel_rcwa[k])))
            print(f"res {res:3d}  level {k}  d={d*1e3:4.0f} nm  "
                  f"T0 FDTD {flux/flux_inc:.3f} / RCWA {rcwa.t_to_T0(t_rcwa[k]):.3f}")
        np.savetxt(out, np.array(rows), delimiter=",", fmt="%.6f",
                   header="level,diameter_nm,T0_fdtd,abs_rel_fdtd,phase_rel_fdtd,"
                          "T0_rcwa,abs_rel_rcwa,phase_rel_rcwa", comments="")
    report(resolutions, t_rcwa, rel_rcwa, diams)


def report(resolutions, t_rcwa, rel_rcwa, diams):
    """Summary table and figure. Only phase *differences between levels* matter
    for the device, so the common phase offset between the two methods is
    removed before comparing."""
    lines = ["# 3D FDTD (Meep) vs RCWA for the 8 unit cells", "",
             "Phase: relative to the empty cell; the common offset between the methods "
             "(circular mean of the differences) is removed, because only the phase "
             "steps between levels affect the vortex plate.", ""]
    results = {}
    for res in resolutions:
        rows = np.loadtxt(RES / f"fdtd_cells_res{res}.csv", delimiter=",", skiprows=1)
        diff = np.angle(np.exp(1j * (rows[:, 4] - rows[:, 7])))
        offset = np.angle(np.mean(np.exp(1j * diff)))
        resid = np.degrees(np.angle(np.exp(1j * (diff - offset))))
        dT = rows[:, 2] - rows[:, 5]
        results[res] = (rows, offset, resid)
        lines += [f"## Resolution {res} px/µm (common offset {np.degrees(offset):+.1f}°)", "",
                  "| Level | d (nm) | T0 FDTD | T0 RCWA | Phase difference after offset (°) |",
                  "|---|---|---|---|---|"]
        lines += [f"| {int(r[0])} | {r[1]:.0f} | {r[2]:.3f} | {r[5]:.3f} | {e:+.1f} |"
                  for r, e in zip(rows, resid)]
        lines += ["", f"RMS phase difference {np.sqrt(np.mean(resid**2)):.1f}°, "
                  f"max |ΔT0| {np.max(np.abs(dT)):.3f}", ""]
    (RES / "fdtd_vs_rcwa.md").write_text("\n".join(lines))
    print("\n".join(lines))

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.4), dpi=150)
    k = np.arange(len(diams))
    ref_phase = np.unwrap(np.angle(rel_rcwa))
    a1.plot(k, ref_phase / (2 * np.pi), "s-", color="#1f4e79", label="RCWA")
    a2.plot(k, rcwa.t_to_T0(t_rcwa), "s-", color="#1f4e79", label="RCWA")
    for (res, (rows, offset, resid)), st in zip(results.items(), ["o--", "^:"]):
        ph = ref_phase + np.radians(resid)
        a1.plot(k, ph / (2 * np.pi), st, color="#7a1530", mfc="none",
                label=f"FDTD {res} px/µm")
        a2.plot(k, rows[:, 2], st, color="#7a1530", mfc="none", label=f"FDTD {res} px/µm")
    a1.set_xlabel("Phase level"); a1.set_ylabel("Phase / 2π (common offset removed)")
    a2.set_xlabel("Phase level"); a2.set_ylabel("Transmittance T0"); a2.set_ylim(0.8, 1.02)
    a1.set_title("Transmitted phase"); a2.set_title("Transmission")
    for a in (a1, a2):
        a.set_xticks(k)
        a.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(RES / "fdtd_vs_rcwa.png")


if __name__ == "__main__":
    main()
