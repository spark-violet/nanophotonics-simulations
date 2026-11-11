"""Transfer-matrix method (TMM) for planar multilayer stacks.

Characteristic-matrix formulation (Born & Wolf / Macleod) for s- and
p-polarised light at arbitrary angle of incidence. Complex refractive
indices n + i*kappa (kappa > 0 for absorption, exp(-i w t) time convention,
the same as Meep) are supported, so absorbing metals such as Ag can be used.

Convention: light comes from the first medium (index n[0]) and leaves into
the last medium (index n[-1]). Both outer media are semi-infinite, so only
the inner layers need a thickness.
"""
from __future__ import annotations

import numpy as np


def _admittance(n: complex, cos_t: complex, pol: str) -> complex:
    """Tilted optical admittance of a layer (in units of free-space admittance)."""
    if pol == "s":
        return n * cos_t
    if pol == "p":
        return n / cos_t
    raise ValueError("pol must be 's' or 'p'")


def amplitudes(n_list, d_list, wavelength, theta0=0.0, pol="s"):
    """Complex amplitude coefficients (r, t) of a multilayer.

    Same inputs as :func:`rt`. The phases are referenced to the first and the
    last interface of the stack.
    """
    return _solve(n_list, d_list, wavelength, theta0, pol)[:2]


def rt(n_list, d_list, wavelength, theta0=0.0, pol="s"):
    """Reflection and transmission of a multilayer at one wavelength and angle.

    Parameters
    ----------
    n_list : sequence of complex
        Refractive indices [incident medium, layer 1, ..., layer N, exit medium].
    d_list : sequence of float
        Thicknesses of the N inner layers, same length unit as ``wavelength``.
    wavelength : float
        Vacuum wavelength.
    theta0 : float
        Angle of incidence in the incident medium (radians).
    pol : {"s", "p"}

    Returns
    -------
    R, T, A : float
        Reflectance, transmittance and absorptance (A = 1 - R - T).
    """
    r, t, eta0, etaN = _solve(n_list, d_list, wavelength, theta0, pol)
    R = abs(r) ** 2
    T = np.real(etaN) / np.real(eta0) * abs(t) ** 2
    return R, T, 1 - R - T


def _solve(n_list, d_list, wavelength, theta0, pol):
    n = np.asarray(n_list, dtype=complex)
    d = np.asarray(d_list, dtype=float)
    if len(d) != len(n) - 2:
        raise ValueError("need one thickness per inner layer")

    # Snell's law with complex angles: n_j sin(theta_j) = n_0 sin(theta_0)
    kx = n[0] * np.sin(theta0)
    cos_t = np.sqrt(1 - (kx / n) ** 2 + 0j)
    # choose the branch with decaying/forward-propagating field
    cos_t = np.where(np.imag(n * cos_t) < 0, -cos_t, cos_t)

    eta = np.array([_admittance(n[j], cos_t[j], pol) for j in range(len(n))])

    M = np.eye(2, dtype=complex)
    for j in range(1, len(n) - 1):
        delta = 2 * np.pi * n[j] * d[j - 1] * cos_t[j] / wavelength
        # time convention exp(-i w t), absorbing media have Im(n) > 0
        layer = np.array(
            [[np.cos(delta), -1j * np.sin(delta) / eta[j]],
             [-1j * eta[j] * np.sin(delta), np.cos(delta)]]
        )
        M = M @ layer

    B, C = M @ np.array([1.0, eta[-1]])
    r = (eta[0] * B - C) / (eta[0] * B + C)
    t = 2 * eta[0] / (eta[0] * B + C)
    return r, t, eta[0], eta[-1]




def spectrum(n_func, d_list, wavelengths, theta0=0.0, pol="s"):
    """R, T, A over a wavelength array. ``n_func(wl)`` returns the index list."""
    out = np.array([rt(n_func(wl), d_list, wl, theta0, pol) for wl in wavelengths])
    return out[:, 0], out[:, 1], out[:, 2]
