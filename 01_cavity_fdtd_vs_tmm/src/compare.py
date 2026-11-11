"""Compare the FDTD (Meep) spectra with the transfer-matrix result.

Outputs (in ../results):
  tmm_normal.csv            TMM R, T, A at normal incidence
  fdtd_vs_tmm.png           overlay of R and T, plus residuals
  convergence.png           RMS difference to TMM vs FDTD grid resolution
  tmm_dispersion.png        TMM angle-resolved reflectance (s and p)
  summary.md                table of resonance wavelength, linewidth and Q
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import stack
import tmm

RES = Path(__file__).resolve().parents[1] / "results"
plt.rcParams.update({"font.size": 10, "figure.dpi": 150})


def resonance(wl, R):
    """Wavelength of the reflectance dip, its FWHM and Q = lambda / FWHM."""
    i = np.argmin(R)
    # refine minimum with a parabola through 3 points
    if 0 < i < len(R) - 1:
        y0, y1, y2 = R[i - 1], R[i], R[i + 1]
        denom = y0 - 2 * y1 + y2
        shift = 0.5 * (y0 - y2) / denom if denom != 0 else 0.0
        lam0 = np.interp(i + shift, np.arange(len(wl)), wl)
    else:
        lam0 = wl[i]
    # FWHM of the dip, measured relative to the local background
    background = 0.5 * (R[: max(1, i // 4)].mean() + R[-max(1, (len(R) - i) // 4):].mean())
    background = max(background, R[i])
    half = R[i] + 0.5 * (background - R[i])
    left = np.where(R[:i] > half)[0]
    right = np.where(R[i:] > half)[0]
    if len(left) == 0 or len(right) == 0:
        return lam0, np.nan, np.nan
    l, r = left[-1], i + right[0]
    wl_l = np.interp(half, [R[l + 1], R[l]], [wl[l + 1], wl[l]])
    wl_r = np.interp(half, [R[r - 1], R[r]], [wl[r - 1], wl[r]])
    fwhm = wl_r - wl_l
    return lam0, fwhm, lam0 / fwhm


def load_fdtd():
    runs = {}
    for f in sorted(RES.glob("fdtd_res*.csv"), key=lambda p: int(p.stem[8:])):
        res = int(f.stem[8:])
        d = np.loadtxt(f, delimiter=",", skiprows=1)
        runs[res] = d
    if not runs:
        raise SystemExit("No FDTD results found. Run run_fdtd.py first.")
    return runs


def main():
    wl = np.linspace(stack.WL_MIN, stack.WL_MAX, 2001)
    R_tmm, T_tmm, A_tmm = tmm.spectrum(stack.index_list, stack.THICKNESSES, wl)
    np.savetxt(RES / "tmm_normal.csv", np.column_stack([wl, R_tmm, T_tmm, A_tmm]),
               delimiter=",", header="wavelength_um,R,T,A", comments="")

    runs = load_fdtd()
    best = max(runs)
    d = runs[best]

    # ---- overlay -------------------------------------------------------
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(6.4, 5.6), sharex=True,
                                  gridspec_kw={"height_ratios": [3, 1]})
    ax.plot(wl * 1e3, R_tmm, color="#1f4e79", lw=2, label="R  (TMM)")
    ax.plot(wl * 1e3, T_tmm, color="#7a1530", lw=2, label="T  (TMM)")
    ax.plot(d[:, 0] * 1e3, d[:, 1], "o", ms=2.5, mfc="none", color="#1f4e79",
            label=f"R  (FDTD, {best} px/µm)")
    ax.plot(d[:, 0] * 1e3, d[:, 2], "o", ms=2.5, mfc="none", color="#7a1530",
            label=f"T  (FDTD, {best} px/µm)")
    ax.set_ylabel("Reflectance / Transmittance")
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, fontsize=8, loc="center right")
    ax.set_title("Ag (30 nm) / spacer n=1.65 (150 nm) / Ag (40 nm) on glass")

    R_i = np.interp(d[:, 0], wl, R_tmm)
    T_i = np.interp(d[:, 0], wl, T_tmm)
    axr.plot(d[:, 0] * 1e3, d[:, 1] - R_i, color="#1f4e79", lw=1, label="ΔR")
    axr.plot(d[:, 0] * 1e3, d[:, 2] - T_i, color="#7a1530", lw=1, label="ΔT")
    axr.axhline(0, color="k", lw=0.5)
    axr.set_xlabel("Wavelength (nm)")
    axr.set_ylabel("FDTD − TMM")
    axr.legend(frameon=False, fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(RES / "fdtd_vs_tmm.png")
    plt.close(fig)

    # ---- resonance metrics and convergence ------------------------------
    lam_t, fwhm_t, q_t = resonance(wl, R_tmm)
    rows = [("TMM", lam_t, fwhm_t, q_t, 0.0, 0.0)]
    res_list, err_list = [], []
    for res, dd in runs.items():
        lam, fwhm, q = resonance(dd[:, 0], dd[:, 1])
        rmse = np.sqrt(np.mean((dd[:, 1] - np.interp(dd[:, 0], wl, R_tmm)) ** 2))
        rows.append((f"FDTD {res} px/µm", lam, fwhm, q, (lam - lam_t) * 1e3, rmse))
        res_list.append(res)
        err_list.append(rmse)

    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    ax.loglog(res_list, err_list, "o-", color="#7a1530", label="FDTD")
    ref = err_list[0] * (np.array(res_list) / res_list[0]) ** -1.0
    ax.loglog(res_list, ref, "--", color="grey", lw=1, label="∝ 1/resolution")
    ax.set_xticks(res_list)
    ax.set_xticklabels([str(r) for r in res_list])
    ax.minorticks_off()
    ax.set_xlabel("FDTD resolution (pixels / µm)")
    ax.set_ylabel("RMS(R_FDTD − R_TMM)")
    ax.set_title("FDTD converges to the TMM result")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(True, which="both", lw=0.3)
    fig.tight_layout()
    fig.savefig(RES / "convergence.png")
    plt.close(fig)

    # ---- TMM angle-resolved reflectance (what a Fourier-space setup measures)
    angles = np.radians(np.linspace(0, 60, 121))
    wl_d = np.linspace(stack.WL_MIN, stack.WL_MAX, 401)
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.6), sharey=True)
    for ax, pol in zip(axes, ["s", "p"]):
        Rmap = np.array([[tmm.rt(stack.index_list(w), stack.THICKNESSES, w, th, pol)[0]
                          for w in wl_d] for th in angles])
        im = ax.pcolormesh(np.degrees(angles), wl_d * 1e3, Rmap.T, shading="auto",
                           cmap="magma", vmin=0, vmax=1)
        ax.set_xlabel("Angle of incidence (°)")
        ax.set_title(f"TMM reflectance, {pol}-pol")
    axes[0].set_ylabel("Wavelength (nm)")
    fig.colorbar(im, ax=axes, label="R")
    fig.savefig(RES / "tmm_dispersion.png", bbox_inches="tight")
    plt.close(fig)

    # ---- summary table --------------------------------------------------
    lines = ["| Method | λ_res (nm) | FWHM (nm) | Q | Δλ vs TMM (nm) | RMS error in R |",
             "|---|---|---|---|---|---|"]
    for name, lam, fwhm, q, dl, rmse in rows:
        lines.append(f"| {name} | {lam*1e3:.2f} | {fwhm*1e3:.2f} | {q:.1f} | "
                     f"{dl:+.2f} | {rmse:.4f} |")
    table = "\n".join(lines)
    (RES / "summary.md").write_text(table + "\n")
    print(table)


if __name__ == "__main__":
    main()
