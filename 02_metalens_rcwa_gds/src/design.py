"""Design parameters for the TiO2 metalens. Lengths in micrometres.

TiO2 pillars on fused silica follow the well-known ALD-based platform for
visible metalenses (Devlin et al., PNAS 2016; Khorasaninejad et al.,
Science 2016). Refractive indices are fixed values at the design wavelength.
"""
WAVELENGTH = 0.633     # HeNe red

N_PILLAR = 2.40        # amorphous / ALD TiO2 at 633 nm (approx.)
N_SUBSTRATE = 1.457    # fused silica
N_SURROUND = 1.0       # air around the pillars

HEIGHT = 0.800         # pillar height (0.6 um gives < 2*pi phase here)
PERIOD = 0.325         # square lattice; < lambda / n_substrate, so no
                       # diffraction into the glass
D_MIN = 0.080          # smallest pillar diameter (fabrication limit)
D_MAX = 0.280          # largest diameter (keeps a >= 45 nm gap)

N_ORDERS = 151         # Fourier orders used for the library

# Lens
LENS_DIAMETER = 40.0   # um
FOCAL_LENGTH = 40.0    # um  -> NA = sin(atan(20/40)) = 0.45
