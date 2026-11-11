"""Definition of the Ag / dielectric / Ag Fabry-Perot microcavity.

All lengths are in micrometres (Meep's default unit, a = 1 um).

Change the numbers here to match a real sample. Both the TMM and the FDTD
script read this file, so the two models always describe the same stack.
"""
from meep.materials import Ag

# Light comes from air, passes the top mirror, the spacer, the bottom mirror
# and leaves into the glass substrate.
N_AIR = 1.0
N_SPACER = 1.65      # ALD Al2O3-like dielectric, treated as non-dispersive
N_GLASS = 1.50       # substrate

T_TOP_AG = 0.030     # 30 nm top mirror
T_SPACER = 0.150     # 150 nm spacer
T_BOTTOM_AG = 0.040  # 40 nm bottom mirror

WL_MIN, WL_MAX = 0.40, 0.90  # wavelength window (um)

# Silver: Meep's built-in Drude-Lorentz fit (Rakic et al., Appl. Opt. 1998).
# The TMM uses the *same* dispersion model, so any difference between the
# two methods comes from the numerics, not from different material data.
SILVER = Ag


def n_silver(wavelength_um: float) -> complex:
    """Complex refractive index of the Meep Ag model at one wavelength."""
    eps = SILVER.epsilon(1.0 / wavelength_um)[0, 0]
    return complex(eps) ** 0.5


def index_list(wavelength_um: float):
    """[incident, top Ag, spacer, bottom Ag, substrate] at one wavelength."""
    n_ag = n_silver(wavelength_um)
    return [N_AIR, n_ag, N_SPACER, n_ag, N_GLASS]


THICKNESSES = [T_TOP_AG, T_SPACER, T_BOTTOM_AG]
