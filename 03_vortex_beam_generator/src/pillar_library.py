"""RCWA pillar library at the geometry chosen by the parameter search, and the
8 pillars used for the 8 phase levels of the vortex plate.

Outputs (in ../results):
  pillar_library.csv     diameter_nm, T0, phase_rad, t_re, t_im
  levels_pillars.csv     (8 chosen pillars, away from sharp resonances) level, target_phase_rad, diameter_nm, T0, phase_rad, t_re, t_im
  pillar_library.png
"""
from pathlib import Path

import grcwa
import matplotlib.pyplot as plt
import numpy as np

import config as C
import rcwa
from parameter_search import best_levels

grcwa.set_backend("numpy")
RES = Path(__file__).resolve().parents[1] / "results"
N_ORDERS = 151
NGRID = 128


def load_geometry():
    g = dict(line.split("=") for line in (RES / "best_geometry.txt").read_text().split())
    return float(g["height_um"]), float(g["period_um"])


def _runs(x, mask):
    """Start/end x of contiguous True runs in mask."""
    out, start = [], None
    for xi, m in zip(x, mask):
        if m and start is None:
            start = xi
        if not m and start is not None:
            out.append((start, prev))
            start = None
        prev = xi
    if start is not None:
        out.append((start, x[-1]))
    return out


def main():
    H, P = load_geometry()
    diams = np.round(np.arange(C.MIN_FEATURE, P - C.MIN_GAP + 1e-9, 0.002), 4)
    t = np.array([rcwa.solve(rcwa.disk_eps(d, P, NGRID), H, P, N_ORDERS)[1] for d in diams])
    T0 = rcwa.t_to_T0(t)
    phase = np.unwrap(np.angle(t))
    np.savetxt(RES / "pillar_library.csv",
               np.column_stack([diams * 1e3, T0, phase, t.real, t.imag]),
               delimiter=",", header="diameter_nm,T0,phase_rad,t_re,t_im", comments="",
               fmt="%.6f")

    # Robustness: skip pillars within +-TOL of a sharp resonance. A diameter is
    # "robust" if T0 stays above 0.9 for every diameter within +-TOL, so a
    # fabrication error of that size cannot land on a resonance.
    TOL = 0.010
    robust = np.array([T0[np.abs(diams - d) <= TOL + 1e-9].min() >= 0.9 for d in diams])
    fom_any, _ = best_levels(t)
    fom_r, idx_r = best_levels(np.where(robust, t, 0))
    idx = idx_r
    fom = fom_r
    print(f"8-level FOM: {fom_any:.3f} using any pillar, {fom_r:.3f} using only "
          f"pillars >= {TOL*1e3:.0f} nm away from resonances")
    # target phases: the chosen pillars define the offset of level 0
    offset = np.angle(np.mean(t[idx] * np.exp(-1j * C.PHASE_LEVELS)))
    targets = C.PHASE_LEVELS + offset
    rows = [(k, targets[k], diams[i] * 1e3, T0[i], np.angle(t[i]), t[i].real, t[i].imag)
            for k, i in enumerate(idx)]
    np.savetxt(RES / "levels_pillars.csv", np.array(rows), delimiter=",", fmt="%.6f",
               header="level,target_phase_rad,diameter_nm,T0,phase_rad,t_re,t_im",
               comments="")
    print(f"H = {H*1e3:.0f} nm, P = {P*1e3:.0f} nm")
    for k, tgt, d, T, ph, *_ in rows:
        err = np.angle(np.exp(1j * (ph - tgt)))
        print(f"  level {k}: d = {d:5.0f} nm   T0 = {T:.3f}   phase error = {np.degrees(err):+5.1f} deg")

    fig, ax = plt.subplots(figsize=(5.8, 3.6), dpi=150)
    ax.plot(diams * 1e3, (phase - phase[0]) / (2 * np.pi), color="#7a1530", lw=2)
    sel_phase = (phase[idx] - phase[0]) / (2 * np.pi)
    ax.plot(diams[idx] * 1e3, sel_phase, "o", mfc="white", mec="#7a1530", ms=6,
            label="8 chosen pillars")
    ax.set_xlabel("Pillar diameter (nm)")
    ax.set_ylabel("Phase / 2π", color="#7a1530")
    ax2 = ax.twinx()
    ax2.plot(diams * 1e3, T0, color="#1f4e79", lw=1.5)
    ax2.set_ylim(0, 1.05)
    ax2.set_ylabel("Transmittance T0", color="#1f4e79")
    for d0, d1 in _runs(diams, ~robust):
        ax.axvspan(d0 * 1e3, d1 * 1e3, color="grey", alpha=0.18, lw=0)
    ax.plot([], [], "s", color="grey", alpha=0.3, label="avoided (near resonance)")
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    ax.set_title(f"TiO$_2$ pillars, H = {H*1e3:.0f} nm, P = {P*1e3:.0f} nm (RCWA)",
                 fontsize=9)
    fig.tight_layout()
    fig.savefig(RES / "pillar_library.png")


if __name__ == "__main__":
    main()
