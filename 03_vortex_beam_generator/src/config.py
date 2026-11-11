"""Shared parameters for the vortex-beam (OAM) metasurface. Lengths in um.

Platform: TiO2 nanostructures on fused silica, illuminated through the
substrate at 633 nm (same platform as project 02).
"""
import numpy as np

WAVELENGTH = 0.633
N_PILLAR = 2.40          # TiO2 (ALD), treated as lossless and non-dispersive
N_SUBSTRATE = 1.457      # fused silica
N_AIR = 1.0

# fabrication limits used in the parameter search and the inverse design
MIN_FEATURE = 0.080      # smallest pillar diameter / smallest structure
MIN_GAP = 0.060          # smallest gap between neighbouring structures
MAX_ASPECT = 10.0        # height / MIN_FEATURE for etching or ALD filling

N_LEVELS = 8             # phase levels of the discretised vortex plate
PHASE_LEVELS = 2 * np.pi * np.arange(N_LEVELS) / N_LEVELS

# vortex lens: phase = l*theta + focusing term
TOPOLOGICAL_CHARGE = 1
PLATE_DIAMETER = 30.0
FOCAL_LENGTH = 30.0
