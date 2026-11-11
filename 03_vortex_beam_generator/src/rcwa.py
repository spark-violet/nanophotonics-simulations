"""Thin wrapper around grcwa for one periodic unit cell.

Light enters from the glass substrate and leaves into air (normal incidence,
x-polarised). The returned t is the complex zeroth-order transmission
coefficient in the usual Fresnel convention, t = n_substrate * t_grcwa
(the same scaling that was validated against the TMM in project 02), so the
power transmittance of the zeroth order is T0 = (n_air / n_substrate) |t|^2.

Works with grcwa's numpy backend and with its autograd backend (used for
the inverse design).
"""
import grcwa
import numpy as np

import config as C


def disk_eps(diameter, period, ngrid):
    """Permittivity grid of a centred circular pillar."""
    x = (np.arange(ngrid) + 0.5) / ngrid * period - period / 2
    X, Y = np.meshgrid(x, x, indexing="ij")
    inside = X**2 + Y**2 <= (diameter / 2) ** 2
    return np.where(inside, C.N_PILLAR**2, C.N_AIR**2).astype(float)


def solve(eps_grid, height, period, n_orders=101):
    """Return (T_total, complex t) for one unit cell.

    eps_grid : 2D array (ngrid x ngrid) of permittivity of the patterned layer
    """
    ngrid = eps_grid.shape[0]
    obj = grcwa.obj(n_orders, [period, 0], [0, period], 1 / C.WAVELENGTH,
                    0.0, 0.0, verbose=0)
    obj.Add_LayerUniform(0.0, C.N_SUBSTRATE**2)
    obj.Add_LayerGrid(height, ngrid, ngrid)
    obj.Add_LayerUniform(0.0, C.N_AIR**2)
    obj.Init_Setup()
    obj.GridLayer_geteps(eps_grid.flatten())
    obj.MakeExcitationPlanewave(1, 0, 0, 0, order=0)
    _, T = obj.RT_Solve(normalize=1)
    E, _ = obj.Solve_FieldFourier(2, 0.0)
    i0 = next(k for k, g in enumerate(obj.G) if g[0] == 0 and g[1] == 0)
    return T, C.N_SUBSTRATE * E[0][i0]


def t_to_T0(t):
    """Zeroth-order power transmittance from the Fresnel-normalised t."""
    return C.N_AIR / C.N_SUBSTRATE * np.abs(t) ** 2
