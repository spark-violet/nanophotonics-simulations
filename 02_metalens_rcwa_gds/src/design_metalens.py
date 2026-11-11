"""Design a metalens from the RCWA library and check how well it focuses.

Steps
  1. Target phase of an ideal lens (hyperbolic profile, focus in air):
         phi(r) = -(2*pi/lambda) * (sqrt(r^2 + f^2) - f)
  2. For each lattice site, choose the pillar whose complex transmission
     t_i is closest to A * exp(i*phi) (A = mean |t|). Only one 2*pi period of
     the library is used, and the complex error (instead of the phase alone)
     also penalises pillars with low transmission.
  3. Build the transmitted field from the chosen pillars (local periodic
     approximation) and propagate it with the angular spectrum method.
  4. Compare the focal spot with an ideal lens of the same aperture.

The propagation is scalar and uses the local periodic approximation, so this
is a fast design check, not a full-wave 3D simulation of the whole lens.

Outputs (in ../results): pillars.csv, metalens_design.png, focusing.png,
metalens_summary.md
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import design

RES = Path(__file__).resolve().parents[1] / "results"
plt.rcParams.update({"font.size": 9, "figure.dpi": 150})


def load_library():
    lib = np.loadtxt(RES / "library.csv", delimiter=",", skiprows=1)
    # keep one 2*pi period of the (unwrapped) phase, so neighbouring pillars
    # change size smoothly; this also excludes the resonance above ~240 nm
    lib = lib[lib[:, 2] < 2 * np.pi]
    d = lib[:, 0] / 1e3
    t = lib[:, 3] + 1j * lib[:, 4]
    t = t * np.exp(-1j * np.angle(t[0]))  # phase relative to thinnest pillar
    return d, t


def target_phase(r):
    k = 2 * np.pi / design.WAVELENGTH
    f = design.FOCAL_LENGTH
    return -k * (np.sqrt(r**2 + f**2) - f)


def lattice():
    P, R = design.PERIOD, design.LENS_DIAMETER / 2
    n = int(np.floor(R / P))
    idx = np.arange(-n, n + 1)
    X, Y = np.meshgrid(idx * P, idx * P, indexing="ij")
    inside = X**2 + Y**2 <= R**2
    return X, Y, inside


def choose_pillars(phi, d_lib, t_lib):
    A = np.mean(np.abs(t_lib))
    target = A * np.exp(1j * phi)
    err = np.abs(t_lib[None, :] - target.reshape(-1, 1))
    k = np.argmin(err, axis=1)
    return d_lib[k], t_lib[k]


def angular_spectrum(u0, dx, wl, z_list):
    """Propagate a 2D field u0 (sampled at dx) to distances z (scalar, air)."""
    n = u0.shape[0]
    fx = np.fft.fftfreq(n, dx)
    FX, FY = np.meshgrid(fx, fx, indexing="ij")
    kz = 2 * np.pi * np.sqrt((1 / wl) ** 2 - FX**2 - FY**2 + 0j)
    kz = np.where(np.imag(kz) < 0, np.conj(kz), kz)   # decaying evanescent
    U0 = np.fft.fft2(u0)
    for z in z_list:
        yield np.fft.ifft2(U0 * np.exp(1j * kz * z))


def fwhm_1d(x, y):
    y = y / y.max()
    i = np.argmax(y)
    left = np.where(y[:i] < 0.5)[0][-1]
    right = i + np.where(y[i:] < 0.5)[0][0]
    xl = np.interp(0.5, [y[left], y[left + 1]], [x[left], x[left + 1]])
    xr = np.interp(0.5, [y[right], y[right - 1]], [x[right], x[right - 1]])
    return xr - xl


def main():
    d_lib, t_lib = load_library()
    X, Y, inside = lattice()
    r = np.hypot(X, Y)
    # A constant phase offset does not change the focus. Pi moves the first
    # 2*pi wrap away from the lens centre, so the central pillars vary smoothly.
    phi = np.mod(target_phase(r) + np.pi, 2 * np.pi)

    d_sel, t_sel = choose_pillars(phi[inside], d_lib, t_lib)
    np.savetxt(RES / "pillars.csv",
               np.column_stack([X[inside], Y[inside], d_sel * 1e3]),
               delimiter=",", header="x_um,y_um,diameter_nm", comments="", fmt="%.4f")

    # ---- transmitted field on a fine grid -------------------------------
    up = 4                                 # sub-samples per unit cell
    dx = design.PERIOD / up
    t_cell = np.zeros(X.shape, complex)
    t_cell[inside] = t_sel
    # ideal lens: perfect phase, same mean amplitude as the chosen pillars,
    # so the Strehl ratio measures the effect of phase and amplitude errors
    ideal_cell = np.where(inside, np.mean(np.abs(t_sel)) * np.exp(1j * target_phase(r)), 0)
    pad = 4 * X.shape[0]                    # zero padding around the lens
    def to_grid(c):
        g = np.kron(c, np.ones((up, up)))
        N = int(2 ** np.ceil(np.log2(g.shape[0] + 2 * pad)))
        out = np.zeros((N, N), complex)
        s = (N - g.shape[0]) // 2
        out[s:s + g.shape[0], s:s + g.shape[1]] = g
        return out
    u_lens, u_ideal = to_grid(t_cell), to_grid(ideal_cell)
    N = u_lens.shape[0]
    x = (np.arange(N) - N / 2 + 0.5) * dx

    f = design.FOCAL_LENGTH
    z_list = np.linspace(0.5 * f, 1.5 * f, 101)
    mid = N // 2
    xz_lens, xz_ideal = [], []
    for ul, ui in zip(angular_spectrum(u_lens, dx, design.WAVELENGTH, z_list),
                      angular_spectrum(u_ideal, dx, design.WAVELENGTH, z_list)):
        xz_lens.append(np.abs(ul[:, mid]) ** 2)
        xz_ideal.append(np.abs(ui[:, mid]) ** 2)
    xz_lens, xz_ideal = np.array(xz_lens), np.array(xz_ideal)
    z_focus = z_list[np.argmax(xz_lens.max(axis=1))]

    u_f = next(angular_spectrum(u_lens, dx, design.WAVELENGTH, [z_focus]))
    u_fi = next(angular_spectrum(u_ideal, dx, design.WAVELENGTH, [f]))
    I_f, I_fi = np.abs(u_f) ** 2, np.abs(u_fi) ** 2
    cut, cut_i = I_f[:, mid], I_fi[:, mid]

    NA = np.sin(np.arctan(design.LENS_DIAMETER / 2 / f))
    fwhm = fwhm_1d(x, cut)
    fwhm_i = fwhm_1d(x, cut_i)
    strehl = I_f.max() / I_fi.max()
    # focusing efficiency: power within a disk of radius 3*FWHM, divided by
    # the power incident on the lens aperture (unit amplitude)
    XX, YY = np.meshgrid(x, x, indexing="ij")
    disk = XX**2 + YY**2 <= (3 * fwhm) ** 2
    # grcwa amplitudes: t_Fresnel = n_sub * t, power T = (n_air / n_sub) |t_Fresnel|^2
    scale = design.N_SURROUND / design.N_SUBSTRATE * design.N_SUBSTRATE**2
    eff = scale * I_f[disk].sum() * dx**2 / (inside.sum() * design.PERIOD**2)

    summary = [
        "# Metalens design summary", "",
        f"- Wavelength: {design.WAVELENGTH*1e3:.0f} nm",
        f"- Diameter: {design.LENS_DIAMETER:.0f} µm, focal length: {f:.0f} µm, NA = {NA:.2f}",
        f"- Pillars: {inside.sum()} TiO2 cylinders, {d_sel.min()*1e3:.0f}–{d_sel.max()*1e3:.0f} nm, "
        f"height {design.HEIGHT*1e3:.0f} nm, period {design.PERIOD*1e3:.0f} nm",
        f"- Focus found at z = {z_focus:.2f} µm (design {f:.0f} µm)",
        f"- Focal spot FWHM: {fwhm*1e3:.0f} nm (ideal lens, same aperture: {fwhm_i*1e3:.0f} nm; "
        f"Airy limit 0.514 λ/NA = {0.514*design.WAVELENGTH/NA*1e3:.0f} nm)",
        f"- Strehl ratio vs ideal lens: {strehl:.3f}",
        f"- Focusing efficiency (power within 3×FWHM / power through aperture): {eff:.2f}",
        "", "Method: RCWA unit-cell library, local periodic approximation, scalar "
        "angular-spectrum propagation.",
    ]
    (RES / "metalens_summary.md").write_text("\n".join(summary) + "\n")
    print("\n".join(summary))

    # ---- figures ---------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.2))
    D = np.full(X.shape, np.nan)
    D[inside] = d_sel * 1e3
    im = axes[0].imshow(D.T, origin="lower", cmap="viridis",
                        extent=[X.min(), X.max(), Y.min(), Y.max()])
    axes[0].set_title("Pillar diameter map (nm)")
    axes[0].set_xlabel("x (µm)"); axes[0].set_ylabel("y (µm)")
    fig.colorbar(im, ax=axes[0], fraction=0.046)

    row = X.shape[0] // 2
    sel = inside[:, row]
    t_row = np.where(sel, t_cell[:, row], np.nan)
    axes[1].plot(X[sel, row], phi[sel, row] / (2 * np.pi), color="grey", lw=1,
                 label="target")
    axes[1].plot(X[sel, row], np.mod(np.angle(t_row[sel]), 2 * np.pi) / (2 * np.pi),
                 ".", ms=2, color="#7a1530", label="implemented")
    axes[1].set_title("Phase profile through the centre")
    axes[1].set_xlabel("x (µm)"); axes[1].set_ylabel("Phase / 2π")
    axes[1].legend(frameon=False, fontsize=7, loc="lower right")

    axes[2].plot(x, cut_i / cut_i.max(), color="grey", lw=1, label="ideal lens")
    axes[2].plot(x, cut / cut_i.max(), color="#7a1530", lw=1.5, label="metalens")
    axes[2].set_xlim(-3, 3)
    axes[2].set_title(f"Focal plane (FWHM {fwhm*1e3:.0f} nm)")
    axes[2].set_xlabel("x (µm)"); axes[2].set_ylabel("Intensity (norm. to ideal)")
    axes[2].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(RES / "metalens_design.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    xr = np.abs(x) <= 8
    ax.imshow(xz_lens[:, xr].T, origin="lower", aspect="auto", cmap="inferno",
              extent=[z_list[0], z_list[-1], x[xr][0], x[xr][-1]])
    ax.axvline(f, color="w", lw=0.6, ls="--")
    ax.set_xlabel("z (µm)"); ax.set_ylabel("x (µm)")
    ax.set_title("Intensity behind the metalens (xz cut)")
    fig.tight_layout()
    fig.savefig(RES / "focusing.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
