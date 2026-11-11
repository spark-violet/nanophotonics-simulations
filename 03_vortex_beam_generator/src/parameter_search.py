"""Parameter search: which pillar height and lattice period work best?

For every (height, period) pair, sweep the pillar diameter with RCWA and ask:
how well can 8 pillars from this sweep realise the 8 phase levels of the
vortex plate, with high transmission?

Figure of merit (FOM), between 0 and 1:

    FOM = | mean_k  sqrt(T0_k) * exp(i (phi_k - target_k)) |^2

where for each level k the pillar is chosen to maximise its contribution and
a common phase offset is optimised. FOM = 1 means 8 perfect, lossless phase
levels. Phase errors and low transmission both reduce it. (FOM is the
zeroth-order efficiency of an 8-level phase plate built from these pillars.)

Fabrication constraints:
  diameters from MIN_FEATURE to period - MIN_GAP;
  height <= MAX_ASPECT * MIN_FEATURE.

The 6 best allowed candidates are then re-evaluated with more Fourier orders
and a finer diameter step (coarse-to-fine search).

Outputs (in ../results): parameter_search.csv, parameter_search_refined.csv,
parameter_search.png, best_geometry.txt
"""
from pathlib import Path

import grcwa
import matplotlib.pyplot as plt
import numpy as np

import config as C
import rcwa

grcwa.set_backend("numpy")
RES = Path(__file__).resolve().parents[1] / "results"
RES.mkdir(exist_ok=True)

HEIGHTS = np.round(np.arange(0.40, 1.001, 0.05), 3)
# upper limit 400 nm: stays below lambda / n_substrate = 434 nm, so the
# reflected light cannot diffract inside the glass
PERIODS = np.round(np.arange(0.250, 0.4001, 0.025), 3)
N_ORDERS = 41      # coarse but fast; the chosen geometry is refined later
NGRID = 64
D_STEP = 0.010


def best_levels(t):
    """FOM and chosen indices for 8 levels from a set of complex t."""
    a = np.sqrt(rcwa.t_to_T0(t)) * np.exp(1j * np.angle(t))
    best = (0.0, None)
    for offset in np.linspace(0, 2 * np.pi, 72, endpoint=False):
        targets = C.PHASE_LEVELS + offset
        proj = np.real(a[None, :] * np.exp(-1j * targets[:, None]))  # level x pillar
        idx = np.argmax(proj, axis=1)
        fom = abs(np.mean(a[idx] * np.exp(-1j * targets))) ** 2
        if fom > best[0]:
            best = (fom, idx)
    return best


def main():
    rows = []
    fom_map = np.full((len(HEIGHTS), len(PERIODS)), np.nan)
    for j, P in enumerate(PERIODS):
        diams = np.arange(C.MIN_FEATURE, P - C.MIN_GAP + 1e-9, D_STEP)
        for i, H in enumerate(HEIGHTS):
            t = np.array([rcwa.solve(rcwa.disk_eps(d, P, NGRID), H, P, N_ORDERS)[1]
                          for d in diams])
            phase = np.unwrap(np.angle(t))
            coverage = (phase.max() - phase.min()) / (2 * np.pi)
            fom, _ = best_levels(t)
            allowed = H <= C.MAX_ASPECT * C.MIN_FEATURE + 1e-9
            fom_map[i, j] = fom
            rows.append((H, P, coverage, np.mean(rcwa.t_to_T0(t)), fom, allowed))
            print(f"H={H:.2f} P={P:.3f}  coverage={coverage:.2f}x2pi  FOM={fom:.3f}"
                  f"{'' if allowed else '  (aspect ratio too high)'}")

    data = np.array(rows, dtype=float)
    np.savetxt(RES / "parameter_search.csv", data, delimiter=",", fmt="%.4f",
               header="height_um,period_um,phase_coverage_2pi,mean_T0,fom,allowed",
               comments="")

    # coarse-to-fine: re-evaluate the 6 best allowed candidates with more
    # Fourier orders and a finer diameter step before choosing
    ok = data[data[:, 5] == 1]
    top = ok[np.argsort(ok[:, 4])[::-1][:6]]
    refined = []
    print("\nrefining the best candidates (101 orders, 5 nm step):")
    for H, P in top[:, :2]:
        diams = np.arange(C.MIN_FEATURE, P - C.MIN_GAP + 1e-9, 0.005)
        t = np.array([rcwa.solve(rcwa.disk_eps(d, P, 128), H, P, 101)[1] for d in diams])
        fom, _ = best_levels(t)
        refined.append((H, P, fom))
        print(f"  H={H:.2f} P={P:.3f}  FOM={fom:.3f}")
    refined = np.array(refined)
    np.savetxt(RES / "parameter_search_refined.csv", refined, delimiter=",", fmt="%.4f",
               header="height_um,period_um,fom_refined", comments="")
    H_best, P_best, fom_best = refined[np.argmax(refined[:, 2])]
    (RES / "best_geometry.txt").write_text(
        f"height_um={H_best:.3f}\nperiod_um={P_best:.3f}\nfom={fom_best:.4f}\n")
    print(f"\nbest allowed geometry: H = {H_best*1e3:.0f} nm, P = {P_best*1e3:.0f} nm, "
          f"FOM = {fom_best:.3f}")

    fig, ax = plt.subplots(figsize=(5.6, 4.2), dpi=150)
    im = ax.pcolormesh(PERIODS * 1e3, HEIGHTS * 1e3, fom_map, cmap="viridis",
                       vmin=0.6, vmax=1, shading="nearest")
    h_max = C.MAX_ASPECT * C.MIN_FEATURE * 1e3
    ax.axhline(h_max, color="w", ls="--", lw=1)
    ax.text(PERIODS[0] * 1e3, h_max + 8, f"aspect-ratio limit ({h_max:.0f} nm)",
            color="w", fontsize=7, va="bottom")
    ax.plot(P_best * 1e3, H_best * 1e3, "*", ms=14, mfc="#e8b4bf", mec="k")
    ax.set_xlabel("Lattice period (nm)")
    ax.set_ylabel("Pillar height (nm)")
    ax.set_title("8-level phase-plate figure of merit (RCWA)")
    fig.colorbar(im, ax=ax, label="FOM")
    fig.tight_layout()
    fig.savefig(RES / "parameter_search.png")


if __name__ == "__main__":
    main()
