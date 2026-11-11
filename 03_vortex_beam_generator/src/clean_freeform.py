"""Fabrication clean-up of the inverse-designed cells, then re-check with RCWA.

The optimised binary patterns can contain one-pixel bumps and notches
(10-30 nm) that e-beam lithography would not reproduce. A morphological
opening and closing with a disk of the minimum feature size removes every
solid feature and gap narrower than that. The cleaned cells are evaluated
again with RCWA (151 orders, 5 nm grid) and the +-5 nm edge-bias test, and
these cleaned results are used for the vortex plate and the GDS.

Inputs : freeform_cells_raw.npz, levels_freeform_raw.csv from inverse_design.py
Outputs: freeform_cells.npz, levels_freeform.csv (cleaned), freeform_cleanup.md
"""
from pathlib import Path

import grcwa
import numpy as np
from scipy.ndimage import binary_closing, binary_opening

import config as C
import rcwa
from inverse_design import _wrap, evaluate, score
from pillar_library import load_geometry
from scipy.ndimage import binary_dilation, binary_erosion

grcwa.set_backend("numpy")
RES = Path(__file__).resolve().parents[1] / "results"


def periodic(op, mask, disk):
    n = mask.shape[0]
    return op(np.tile(mask, (3, 3)), structure=disk)[n:2 * n, n:2 * n]


def main():
    H, P = load_geometry()
    z = np.load(RES / "freeform_cells_raw.npz")
    designs, pixel = z["designs"].astype(bool), float(z["pixel"])
    old = np.loadtxt(RES / "levels_freeform_raw.csv", delimiter=",", skiprows=1)
    psis = old[:, 1]

    r = int(np.ceil(C.MIN_FEATURE / pixel / 2))          # disk radius in pixels
    y, x = np.ogrid[-r:r + 1, -r:r + 1]
    disk = x**2 + y**2 <= r**2

    lines = ["# Fabrication clean-up of the freeform cells", "",
             f"Opening + closing with a {(2*r+1)*pixel*1e3:.0f} nm disk removes features "
             f"and gaps below the {C.MIN_FEATURE*1e3:.0f} nm limit. Scores re-computed with RCWA.",
             "", "| Level | Pixels changed | Score before | Score after | T0 after | "
             "Phase error after (°) | Worst ±5 nm after |", "|---|---|---|---|---|---|---|"]
    cleaned, rows = [], []
    for k, rho in enumerate(designs):
        c = periodic(binary_closing, periodic(binary_opening, rho, disk), disk)
        changed = int(np.sum(c != rho))
        fine = np.kron(c, np.ones((2, 2))).astype(bool)
        t = evaluate(fine.astype(float), H, P)
        t_er = evaluate(_wrap(binary_erosion, fine), H, P)
        t_di = evaluate(_wrap(binary_dilation, fine), H, P)
        worst = min(score(t_er, psis[k]), score(t_di, psis[k]))
        err = np.degrees(np.angle(t * np.exp(-1j * psis[k])))
        cleaned.append(c.astype(float))
        rows.append((k, psis[k], rcwa.t_to_T0(t), np.angle(t), t.real, t.imag,
                     score(t, psis[k]), old[k, 7], worst, old[k, 9]))
        lines.append(f"| {k} | {changed} | {old[k, 6]:.3f} | {score(t, psis[k]):.3f} | "
                     f"{rcwa.t_to_T0(t):.3f} | {err:+.1f} | {worst:.3f} |")
        print(lines[-1])

    np.savez(RES / "freeform_cells.npz", designs=np.array(cleaned), period=P, height=H,
             pixel=pixel)
    rows = np.array(rows)
    np.savetxt(RES / "levels_freeform.csv", rows, delimiter=",", fmt="%.6f",
               header="level,target_phase_rad,T0,phase_rad,t_re,t_im,score_freeform,"
                      "score_pillar,worst_freeform_pm5nm,worst_pillar_pm5nm", comments="")

    def fom(t):
        return abs(np.mean(np.sqrt(C.N_AIR / C.N_SUBSTRATE) * t * np.exp(-1j * psis))) ** 2
    lp = np.loadtxt(RES / "levels_pillars.csv", delimiter=",", skiprows=1)
    lines += ["", f"8-level FOM: pillars {fom(lp[:, 5] + 1j * lp[:, 6]):.3f}, "
              f"freeform (cleaned) {fom(rows[:, 4] + 1j * rows[:, 5]):.3f}",
              f"Mean worst-case score (±5 nm): pillars {rows[:, 9].mean():.3f}, "
              f"freeform (cleaned) {rows[:, 8].mean():.3f}"]
    (RES / "freeform_cleanup.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[-2:]))


if __name__ == "__main__":
    main()
