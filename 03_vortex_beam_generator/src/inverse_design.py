"""Inverse design of freeform unit cells for the 8 phase levels.

Instead of choosing from a family of circular pillars, the shape inside each
unit cell is optimised directly. The cell keeps the same height and period
as the pillars, so it could be made in the same fabrication step.

Method (density-based topology optimisation):
  - design grid of 40 x 40 pixels (10 nm) inside the 400 nm cell
  - C4v symmetry (mirror x, mirror y, diagonal), so the cell responds the
    same to x- and y-polarised light at normal incidence
  - density filter (Gaussian, sigma 40 nm, periodic) to suppress small features,
    then a tanh projection whose steepness beta is increased in stages so
    the design ends up binary (TiO2 or air)
  - objective for level k:  J = sqrt(n_air / n_sub) * Re( t * exp(-i psi_k) )
                              = sqrt(T0) * cos(phase error)
    which is 1 for a lossless cell with exactly the target phase
  - gradients of J with respect to all pixels from automatic differentiation
    through the RCWA solver (grcwa with the autograd backend), optimiser
    L-BFGS-B with bounds [0, 1]
  - robust formulation: the objective is a soft minimum over the nominal
    design and an eroded and a dilated version (thresholds 0.35/0.5/0.65)
  - each level is started from its chosen pillar (levels_pillars.csv)

The final binary designs are re-evaluated with the numpy backend at higher
accuracy (151 Fourier orders, 5 nm grid), and checked for robustness by
moving every edge 5 nm in and out. The same check is applied to the pillars.

Outputs (in ../results): freeform_cells_raw.npz, levels_freeform_raw.csv,
inverse_design.png, inverse_design_log.md   (clean_freeform.py then makes the
fabrication-cleaned versions used downstream)
"""
from pathlib import Path

import autograd.numpy as anp
import grcwa
import matplotlib.pyplot as plt
import numpy as np
from autograd import value_and_grad
from scipy.ndimage import binary_dilation, binary_erosion
from scipy.optimize import minimize

import config as C
import rcwa
from pillar_library import load_geometry

RES = Path(__file__).resolve().parents[1] / "results"
NPIX = 40
FILTER_SIGMA_PX = 4.0          # Gaussian filter width (pixels of 10 nm)
BETAS = [4, 8, 16, 32]
ITERS_PER_BETA = 30
N_ORDERS_OPT = 101
N_ORDERS_CHECK = 151
# robust optimisation: the projection threshold eta = 0.5 -/+ 0.15 gives a
# dilated / eroded version of the same design (edges shifted by ~10 nm).
# The objective is a soft minimum over the three, so the optimiser improves
# the worst of them and cannot rely on a sharp resonance that only works
# for one exact shape. Together with the filter, the eroded and dilated
# designs also impose a minimum length scale on solid and void.
ETAS = [0.35, 0.5, 0.65]
SOFTMIN_K = 25.0               # sharpness of the soft minimum over the three
EDGE_SHIFT_UPSAMPLE = 2        # check grid: 5 nm pixels


# ---------------------------------------------------------------- helpers --
def symmetrize(x):
    """C4v symmetric field from an arbitrary square array."""
    x = 0.5 * (x + x.T)
    x = 0.5 * (x + x[::-1, :])
    x = 0.5 * (x + x[:, ::-1])
    return x


def _gauss_kernel_fft(n, sigma):
    k = np.fft.fftfreq(n)
    KX, KY = np.meshgrid(k, k, indexing="ij")
    return np.exp(-2 * (np.pi * sigma) ** 2 * (KX**2 + KY**2))


KERNEL = _gauss_kernel_fft(NPIX, FILTER_SIGMA_PX)


def density(x, beta, eta=0.5):
    """design variables -> symmetric, filtered, projected density in [0, 1]."""
    x = symmetrize(x)
    xf = anp.real(anp.fft.ifft2(anp.fft.fft2(x) * KERNEL))   # periodic filter
    num = anp.tanh(beta * eta) + anp.tanh(beta * (xf - eta))
    den = anp.tanh(beta * eta) + anp.tanh(beta * (1 - eta))
    return num / den


def eps_from_density(rho):
    return C.N_AIR**2 + rho * (C.N_PILLAR**2 - C.N_AIR**2)


def score(t, psi):
    return np.sqrt(C.N_AIR / C.N_SUBSTRATE) * np.real(t * np.exp(-1j * psi))


def evaluate(rho_binary, H, P, n_orders=N_ORDERS_CHECK, upsample=1):
    grcwa.set_backend("numpy")
    eps = eps_from_density(np.kron(rho_binary, np.ones((upsample, upsample))))
    return rcwa.solve(eps, H, P, n_orders)[1]


def disk_density(diameter, P):
    x = (np.arange(NPIX) + 0.5) / NPIX * P - P / 2
    X, Y = np.meshgrid(x, x, indexing="ij")
    return (X**2 + Y**2 <= (diameter / 2) ** 2).astype(float)


def _wrap(op, mask):
    """Erosion/dilation with periodic boundaries (the cell repeats)."""
    n = mask.shape[0]
    big = np.tile(mask, (3, 3))
    return op(big, iterations=1)[n:2 * n, n:2 * n].astype(float)


# ------------------------------------------------------------ optimise one --
def optimise_level(psi, x0, H, P):
    grcwa.set_backend("autograd")

    def objective(xflat, beta):
        acc = 0.0
        for eta in ETAS:
            rho = density(xflat.reshape(NPIX, NPIX), beta, eta)
            _, t = rcwa.solve(eps_from_density(rho), H, P, N_ORDERS_OPT)
            Ji = anp.sqrt(C.N_AIR / C.N_SUBSTRATE) * anp.real(t * anp.exp(-1j * psi))
            acc = acc + anp.exp(-SOFTMIN_K * Ji)
        soft_min = -anp.log(acc / len(ETAS)) / SOFTMIN_K
        return -soft_min

    x = x0.flatten()
    history = []
    for beta in BETAS:
        fg = value_and_grad(lambda v: objective(v, beta))
        res = minimize(lambda v: tuple(np.asarray(a, float) for a in fg(v)), x,
                       jac=True, method="L-BFGS-B", bounds=[(0, 1)] * x.size,
                       options={"maxiter": ITERS_PER_BETA})
        x = res.x
        history.append((beta, -res.fun))
    rho = density(x.reshape(NPIX, NPIX), BETAS[-1])
    return (np.asarray(rho) > 0.5).astype(float), history


def gradient_check(H, P):
    """Compare the autograd gradient with a finite difference (one pixel)."""
    grcwa.set_backend("autograd")
    rng = np.random.default_rng(0)
    x = rng.uniform(0.3, 0.7, (NPIX, NPIX))

    def f(v):
        _, t = rcwa.solve(eps_from_density(density(v.reshape(NPIX, NPIX), 4)), H, P, 41)
        return anp.real(t)

    val, g = value_and_grad(f)(x.flatten())
    i = NPIX * 12 + 17
    h = 1e-5
    xp, xm = x.flatten().copy(), x.flatten().copy()
    xp[i] += h
    xm[i] -= h
    fd = (f(xp) - f(xm)) / (2 * h)
    return float(g[i]), float(fd)


# ------------------------------------------------------------------- main --
def main():
    H, P = load_geometry()
    levels = np.loadtxt(RES / "levels_pillars.csv", delimiter=",", skiprows=1)
    psis, d_pillar = levels[:, 1], levels[:, 2] / 1e3
    t_pillar = levels[:, 5] + 1j * levels[:, 6]

    g_ad, g_fd = gradient_check(H, P)
    log = ["# Inverse design log", "",
           f"Gradient check (one pixel): autograd {g_ad:.6e}, finite difference {g_fd:.6e}, "
           f"relative difference {abs(g_ad - g_fd) / abs(g_fd):.1e}", ""]
    print(log[2])

    designs, rows = [], []
    log += ["Score = sqrt(T0)·cos(phase error), 1 is perfect. 'Worst ±5 nm' = lowest score "
            "when every edge moves 5 nm in or out (pillar diameter ±10 nm).", "",
            "| Level | Pillar score | Pillar worst ±5 nm | Freeform score | Freeform T0 | "
            "Freeform phase error (°) | Freeform worst ±5 nm |",
            "|---|---|---|---|---|---|---|"]
    for k, psi in enumerate(psis):
        x0 = 0.15 + 0.7 * disk_density(d_pillar[k], P)   # start from the pillar
        rho, hist = optimise_level(psi, x0, H, P)
        # check on a 5 nm grid with more Fourier orders than in the optimisation
        u = EDGE_SHIFT_UPSAMPLE
        fine = np.kron(rho, np.ones((u, u))).astype(bool)
        t = evaluate(fine.astype(float), H, P)
        # robustness: shift every edge by one 5 nm pixel (erosion and dilation),
        # i.e. a pillar diameter changes by 10 nm
        t_er = evaluate(_wrap(binary_erosion, fine), H, P)
        t_di = evaluate(_wrap(binary_dilation, fine), H, P)
        worst = min(score(t_er, psi), score(t_di, psi))
        p_er = rcwa.solve(rcwa.disk_eps(d_pillar[k] - 0.010, P, 2 * NPIX), H, P, N_ORDERS_CHECK)[1]
        p_di = rcwa.solve(rcwa.disk_eps(d_pillar[k] + 0.010, P, 2 * NPIX), H, P, N_ORDERS_CHECK)[1]
        worst_p = min(score(p_er, psi), score(p_di, psi))
        err = np.degrees(np.angle(t * np.exp(-1j * psi)))
        designs.append(rho)
        rows.append((k, psi, rcwa.t_to_T0(t), np.angle(t), t.real, t.imag,
                     score(t, psi), score(t_pillar[k], psi), worst, worst_p))
        log.append(f"| {k} | {score(t_pillar[k], psi):.3f} | {worst_p:.3f} | {score(t, psi):.3f} | "
                   f"{rcwa.t_to_T0(t):.3f} | {err:+.1f} | {worst:.3f} |")
        print(log[-1], "  beta stages:", " ".join(f"{b}:{j:.3f}" for b, j in hist))

    rows = np.array(rows)
    np.savez(RES / "freeform_cells_raw.npz", designs=np.array(designs), period=P, height=H,
             pixel=P / NPIX)
    np.savetxt(RES / "levels_freeform_raw.csv", rows, delimiter=",", fmt="%.6f",
               header="level,target_phase_rad,T0,phase_rad,t_re,t_im,score_freeform,"
                      "score_pillar,worst_freeform_pm5nm,worst_pillar_pm5nm", comments="")

    def plate_fom(t_levels):
        a = np.sqrt(C.N_AIR / C.N_SUBSTRATE) * t_levels * np.exp(-1j * psis)
        return abs(np.mean(a)) ** 2
    t_free = rows[:, 4] + 1j * rows[:, 5]
    log += ["", f"8-level FOM (zeroth-order efficiency of the phase plate): "
            f"pillars {plate_fom(t_pillar):.3f}, freeform {plate_fom(t_free):.3f}"]
    (RES / "inverse_design_log.md").write_text("\n".join(log) + "\n")
    print(log[-1])

    # ---- figure: the 8 freeform cells and their pillar counterparts --------
    fig, axes = plt.subplots(2, 8, figsize=(12, 3.4), dpi=150)
    for k in range(8):
        axes[0, k].imshow(disk_density(d_pillar[k], P).T, cmap="Greys", origin="lower",
                          extent=[0, P * 1e3, 0, P * 1e3], vmin=0, vmax=1)
        axes[1, k].imshow(designs[k].T, cmap="Greys", origin="lower",
                          extent=[0, P * 1e3, 0, P * 1e3], vmin=0, vmax=1)
        axes[0, k].set_title(f"level {k}\n{rows[k, 7]:.2f}", fontsize=7)
        axes[1, k].set_title(f"{rows[k, 6]:.2f}", fontsize=7)
        for a in axes[:, k]:
            a.set_xticks([]); a.set_yticks([])
    axes[0, 0].set_ylabel("pillar", fontsize=8)
    axes[1, 0].set_ylabel("freeform", fontsize=8)
    fig.suptitle("Unit cells (TiO$_2$ = black, 400 nm × 400 nm). "
                 "Numbers: sqrt(T0)·cos(phase error), 1 = perfect", fontsize=9)
    fig.tight_layout()
    fig.savefig(RES / "inverse_design.png")


if __name__ == "__main__":
    main()
