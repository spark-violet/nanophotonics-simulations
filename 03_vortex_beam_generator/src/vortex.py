"""Build the vortex-beam generator and measure the beam it makes.

The metasurface imposes   phi(x, y) = l * atan2(y, x) - k (sqrt(r^2 + f^2) - f)
i.e. a spiral phase (orbital angular momentum l) combined with a lens, so
the vortex (a ring of light with a phase singularity in the middle) forms
at the focal plane. The phase is discretised into the 8 levels, and each
lattice site gets the unit cell of its level.

Unit-cell sets compared (all use the local periodic approximation):
  ideal_continuous : perfect, lossless, continuous phase
  ideal_8_level    : perfect, lossless, 8 phase levels
  pillars_rcwa     : 8 TiO2 pillars, t from RCWA
  pillars_fdtd     : the same pillars, t from 3D FDTD (if available)
  freeform_inverse : inverse-designed cells, t from RCWA

Beam metrics at the focal plane:
  OAM purity  : power in the target mode l / total power, from an azimuthal
                Fourier decomposition of the field on rings
  efficiency  : power inside a 3 um radius / power incident on the plate

Outputs (in ../results): vortex_metrics.md, vortex_beams.png, levels_used.csv
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import map_coordinates

import config as C
from pillar_library import load_geometry

RES = Path(__file__).resolve().parents[1] / "results"
UP = 4                  # field samples per unit cell
R_WINDOW = 3.0          # um, radius used for efficiency and OAM analysis
M_MAX = 10              # OAM modes analysed: -M_MAX ... M_MAX


# ------------------------------------------------------------ level sets ----
def level_sets():
    lp = np.loadtxt(RES / "levels_pillars.csv", delimiter=",", skiprows=1)
    psi = lp[:, 1]
    sets = {"ideal_8_level": np.sqrt(C.N_SUBSTRATE / C.N_AIR) * np.exp(1j * psi),
            "pillars_rcwa": lp[:, 5] + 1j * lp[:, 6]}
    fdtd_files = sorted(RES.glob("fdtd_cells_res*.csv"), key=lambda p: int(p.stem[14:]))
    if fdtd_files:
        f = np.loadtxt(fdtd_files[-1], delimiter=",", skiprows=1)
        # FDTD gives T0 and the phase relative to the empty cell; use the RCWA
        # empty-cell phase so that both pillar sets share the same reference
        ref = np.angle(sets["pillars_rcwa"]) - f[:, 7]      # RCWA phase of empty cell
        amp = np.sqrt(f[:, 2] * C.N_SUBSTRATE / C.N_AIR)
        sets[f"pillars_fdtd"] = amp * np.exp(1j * (f[:, 4] + ref))
    ff = RES / "levels_freeform.csv"
    if ff.exists():
        f = np.loadtxt(ff, delimiter=",", skiprows=1)
        sets["freeform_inverse"] = f[:, 4] + 1j * f[:, 5]
    return psi, sets


# --------------------------------------------------------- plate and beam ---
def plate_phase(X, Y, charge):
    k = 2 * np.pi / C.WAVELENGTH
    f = C.FOCAL_LENGTH
    return charge * np.arctan2(Y, X) - k * (np.sqrt(X**2 + Y**2 + f**2) - f)


def build_field(P, charge, psi, t_levels, continuous=False):
    n = int(np.floor(C.PLATE_DIAMETER / 2 / P))
    idx = np.arange(-n, n + 1) * P
    X, Y = np.meshgrid(idx, idx, indexing="ij")
    inside = X**2 + Y**2 <= (C.PLATE_DIAMETER / 2) ** 2
    phi = plate_phase(X, Y, charge)
    if continuous:
        cell = np.sqrt(C.N_SUBSTRATE / C.N_AIR) * np.exp(1j * (phi + psi[0]))
    else:
        level = np.mod(np.round(phi / (2 * np.pi / C.N_LEVELS)), C.N_LEVELS).astype(int)
        cell = t_levels[level]
    cell = np.where(inside, cell, 0)
    # scale so |t|^2 * n_air / n_sub is the power transmittance
    cell = cell * np.sqrt(C.N_AIR / C.N_SUBSTRATE)
    g = np.kron(cell, np.ones((UP, UP)))
    N = int(2 ** np.ceil(np.log2(3 * g.shape[0])))
    u = np.zeros((N, N), complex)
    s = (N - g.shape[0]) // 2
    u[s:s + g.shape[0], s:s + g.shape[1]] = g
    return u, P / UP, inside.sum() * P**2


def propagate(u, dx, z):
    n = u.shape[0]
    fx = np.fft.fftfreq(n, dx)
    FX, FY = np.meshgrid(fx, fx, indexing="ij")
    kz = 2 * np.pi * np.sqrt((1 / C.WAVELENGTH) ** 2 - FX**2 - FY**2 + 0j)
    return np.fft.ifft2(np.fft.fft2(u) * np.exp(1j * kz * z))


def oam_spectrum(E, dx):
    """Power in each OAM mode m, integrated over rings r < R_WINDOW."""
    n = E.shape[0]
    c = n / 2 - 0.5
    radii = np.arange(dx, R_WINDOW, dx / 2)
    theta = np.linspace(0, 2 * np.pi, 256, endpoint=False)
    power = np.zeros(2 * M_MAX + 1)
    for r in radii:
        xs = c + r * np.cos(theta) / dx
        ys = c + r * np.sin(theta) / dx
        ring = (map_coordinates(E.real, [xs, ys], order=3)
                + 1j * map_coordinates(E.imag, [xs, ys], order=3))
        cm = np.fft.fft(ring) / len(theta)
        m = np.fft.fftfreq(len(theta), 1 / len(theta)).astype(int)
        for j, mm in enumerate(range(-M_MAX, M_MAX + 1)):
            power[j] += abs(cm[m == mm][0]) ** 2 * r
    return np.arange(-M_MAX, M_MAX + 1), power / power.sum()


def analyse(P, charge, psi, t_levels, continuous=False):
    u, dx, area = build_field(P, charge, psi, t_levels, continuous)
    E = propagate(u, dx, C.FOCAL_LENGTH)
    n = E.shape[0]
    x = (np.arange(n) - n / 2 + 0.5) * dx
    XX, YY = np.meshgrid(x, x, indexing="ij")
    disk = XX**2 + YY**2 <= R_WINDOW**2
    I = np.abs(E) ** 2          # |cell|^2 = power transmittance (see build_field)
    eff = I[disk].sum() * dx**2 / area
    m, pw = oam_spectrum(E, dx)
    purity = pw[m == charge][0]
    return dict(E=E, x=x, eff=eff, purity=purity, m=m, oam=pw)


# --------------------------------------------------- fabrication tolerance ---
def biased_sets(H, P, psi, bias_nm):
    """Unit-cell t when every edge moves by bias_nm (systematic fab bias)."""
    import grcwa
    from scipy.ndimage import binary_dilation, binary_erosion
    import rcwa
    grcwa.set_backend("numpy")
    out = {}
    lp = np.loadtxt(RES / "levels_pillars.csv", delimiter=",", skiprows=1)
    d = lp[:, 2] / 1e3 + 2 * bias_nm / 1e3
    out["pillars_rcwa"] = np.array([rcwa.solve(rcwa.disk_eps(di, P, 80), H, P, 151)[1]
                                    for di in d])
    npz = RES / "freeform_cells.npz"
    if npz.exists():
        z = np.load(npz)
        designs, pix = z["designs"], float(z["pixel"])
        up = int(round(pix / 0.005))            # 5 nm grid
        steps = int(round(abs(bias_nm) / 5))
        op = binary_dilation if bias_nm > 0 else binary_erosion
        ts = []
        for rho in designs:
            fine = np.kron(rho, np.ones((up, up))).astype(bool)
            n = fine.shape[0]
            if steps:
                fine = op(np.tile(fine, (3, 3)), iterations=steps)[n:2 * n, n:2 * n]
            eps = C.N_AIR**2 + fine * (C.N_PILLAR**2 - C.N_AIR**2)
            ts.append(rcwa.solve(eps.astype(float), H, P, 151)[1])
        out["freeform_inverse"] = np.array(ts)
    return out


# -------------------------------------------------------------------- main --
def main():
    H, P = load_geometry()
    psi, sets = level_sets()
    np.savetxt(RES / "levels_used.csv",
               np.column_stack([np.arange(8), psi] + [np.angle(v) for v in sets.values()]
                               + [np.abs(v) for v in sets.values()]),
               delimiter=",", fmt="%.5f",
               header="level,target_phase," + ",".join(f"phase_{k}" for k in sets)
                      + "," + ",".join(f"abs_{k}" for k in sets), comments="")

    names = ["ideal_continuous"] + list(sets)
    lines = ["# Vortex-beam generator: results", "",
             f"Plate: {C.PLATE_DIAMETER:.0f} µm diameter, focal length {C.FOCAL_LENGTH:.0f} µm, "
             f"{C.N_LEVELS} phase levels, TiO2 cells H = {H*1e3:.0f} nm, P = {P*1e3:.0f} nm, "
             f"λ = {C.WAVELENGTH*1e3:.0f} nm.", "",
             "| Unit cells | l | OAM purity | Efficiency (r < 3 µm) |", "|---|---|---|---|"]
    results = {}
    for charge in (1, 2):
        for name in names:
            cont = name == "ideal_continuous"
            t = None if cont else sets[name]
            r = analyse(P, charge, psi, t, continuous=cont)
            results[(name, charge)] = r
            lines.append(f"| {name} | {charge} | {r['purity']:.3f} | {r['eff']:.3f} |")
            print(lines[-1])
    # fabrication tolerance: all edges moved by -5 / +5 nm
    lines += ["", "## Fabrication tolerance (l = 1)", "",
              "Every edge moved by the same amount (systematic over- or under-etch).", "",
              "| Unit cells | Edge bias | OAM purity | Efficiency (r < 3 µm) |",
              "|---|---|---|---|"]
    for bias in (-5, 5):
        for name, t in biased_sets(H, P, psi, bias).items():
            r = analyse(P, 1, psi, t)
            lines.append(f"| {name} | {bias:+d} nm | {r['purity']:.3f} | {r['eff']:.3f} |")
            print(lines[-1])
    lines += ["", "OAM purity: share of the focal-plane power in the target mode "
              "(azimuthal decomposition on rings, r < 3 µm). Efficiency: power within "
              "r < 3 µm divided by the power incident on the plate."]
    (RES / "vortex_metrics.md").write_text("\n".join(lines) + "\n")

    # ---- figure ----------------------------------------------------------
    show = [n for n in ["ideal_continuous", "pillars_rcwa", "freeform_inverse"] if
            (n, 1) in results]
    fig, axes = plt.subplots(len(show), 4, figsize=(12, 2.9 * len(show)), dpi=130)
    axes = np.atleast_2d(axes)
    for row, name in enumerate(show):
        for col, charge in enumerate((1, 2)):
            r = results[(name, charge)]
            sel = np.abs(r["x"]) <= 3
            E = r["E"][np.ix_(sel, sel)]
            ext = [r["x"][sel][0], r["x"][sel][-1]] * 2
            ax = axes[row, 2 * col]
            ax.imshow((np.abs(E) ** 2).T, origin="lower", cmap="inferno", extent=ext)
            ax.set_title(f"{name}, l={charge}\nintensity, purity {r['purity']:.2f}",
                         fontsize=8)
            ax2 = axes[row, 2 * col + 1]
            ax2.imshow(np.angle(E).T, origin="lower", cmap="twilight", extent=ext,
                       vmin=-np.pi, vmax=np.pi)
            ax2.set_title("phase", fontsize=8)
            for a in (ax, ax2):
                a.set_xlabel("x (µm)", fontsize=7)
                a.tick_params(labelsize=6)
    fig.tight_layout()
    fig.savefig(RES / "vortex_beams.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 3.2), dpi=150)
    w = 0.8 / len(names)
    for i, name in enumerate(names):
        r = results[(name, 1)]
        ax.bar(r["m"] + (i - len(names) / 2 + 0.5) * w, r["oam"], width=w, label=name)
    ax.set_xlim(-4.5, 6.5)
    ax.set_yscale("log"); ax.set_ylim(1e-5, 1.5)
    ax.set_xlabel("OAM mode m"); ax.set_ylabel("Power fraction")
    ax.set_title("OAM spectrum at the focus, target l = 1")
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(RES / "oam_spectrum.png")


if __name__ == "__main__":
    main()
