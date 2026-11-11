"""1D FDTD simulation (Meep) of the Ag/dielectric/Ag microcavity at normal incidence.

Reflectance and transmittance are obtained with the standard two-run flux
method: a reference run without the stack records the incident fields, which
are then subtracted from the reflected-flux monitor in the run with the stack.

Usage:  python run_fdtd.py [resolution ...]
        (default resolutions: 200 400 800 pixels/um)
Writes results/fdtd_res<resolution>.csv with columns wavelength_um, R, T.
"""
import sys
from pathlib import Path

import meep as mp
import numpy as np

import stack

mp.verbosity(0)
OUT = Path(__file__).resolve().parents[1] / "results"
OUT.mkdir(exist_ok=True)

# Frequency window (Meep frequency = 1 / wavelength in um)
fmin, fmax = 1 / stack.WL_MAX, 1 / stack.WL_MIN
fcen, df = 0.5 * (fmin + fmax), (fmax - fmin) * 1.2
NFREQ = 400

PML = 1.0
AIR_GAP = 1.0      # space between source and stack
GLASS_GAP = 1.0    # glass below the stack before the PML


def geometry_and_cell():
    t1, ts, t2 = stack.THICKNESSES
    stack_len = t1 + ts + t2
    sz = 2 * PML + AIR_GAP + stack_len + GLASS_GAP
    z_top = 0.5 * sz - PML - AIR_GAP          # top surface of the stack
    z = z_top
    geom = []
    for thick, mat in [(t1, stack.SILVER),
                       (ts, mp.Medium(index=stack.N_SPACER)),
                       (t2, stack.SILVER)]:
        geom.append(mp.Block(size=mp.Vector3(mp.inf, mp.inf, thick),
                             center=mp.Vector3(z=z - thick / 2), material=mat))
        z -= thick
    # glass substrate from the bottom of the stack into the lower PML
    glass_len = z + 0.5 * sz
    geom.append(mp.Block(size=mp.Vector3(mp.inf, mp.inf, glass_len),
                         center=mp.Vector3(z=z - glass_len / 2),
                         material=mp.Medium(index=stack.N_GLASS)))
    return geom, sz, z_top, z


def run(resolution: int):
    geom, sz, z_top, z_bottom = geometry_and_cell()
    cell = mp.Vector3(0, 0, sz)
    pml = [mp.PML(PML)]
    z_src = 0.5 * sz - PML - 0.2
    z_refl = z_src - 0.1
    z_tran = z_bottom - 0.3

    src = [mp.Source(mp.GaussianSource(fcen, fwidth=df),
                     component=mp.Ex, center=mp.Vector3(z=z_src))]

    def sim_with(geometry, default):
        return mp.Simulation(cell_size=cell, boundary_layers=pml, sources=src,
                             geometry=geometry, default_material=default,
                             resolution=resolution, dimensions=1)

    # --- reference run: empty air, transmission monitor also in air
    ref = sim_with([], mp.air)
    refl_ref = ref.add_flux(fcen, df, NFREQ, mp.FluxRegion(center=mp.Vector3(z=z_refl)))
    tran_ref = ref.add_flux(fcen, df, NFREQ, mp.FluxRegion(center=mp.Vector3(z=z_tran)))
    ref.run(until_after_sources=mp.stop_when_fields_decayed(
        20, mp.Ex, mp.Vector3(z=z_tran), 1e-9))
    incident = np.array(mp.get_fluxes(tran_ref))
    refl_data = ref.get_flux_data(refl_ref)
    freqs = np.array(mp.get_flux_freqs(refl_ref))

    # --- run with the cavity
    sim = sim_with(geom, mp.air)
    refl = sim.add_flux(fcen, df, NFREQ, mp.FluxRegion(center=mp.Vector3(z=z_refl)))
    tran = sim.add_flux(fcen, df, NFREQ, mp.FluxRegion(center=mp.Vector3(z=z_tran)))
    sim.load_minus_flux_data(refl, refl_data)
    sim.run(until_after_sources=mp.stop_when_fields_decayed(
        20, mp.Ex, mp.Vector3(z=z_refl), 1e-9))

    R = -np.array(mp.get_fluxes(refl)) / incident
    T = np.array(mp.get_fluxes(tran)) / incident

    wl = 1 / freqs
    keep = (wl >= stack.WL_MIN) & (wl <= stack.WL_MAX)
    order = np.argsort(wl[keep])
    data = np.column_stack([wl[keep][order], R[keep][order], T[keep][order]])
    path = OUT / f"fdtd_res{resolution}.csv"
    np.savetxt(path, data, delimiter=",", header="wavelength_um,R,T", comments="")
    print(f"resolution {resolution:4d} px/um  ->  {path.name}")


if __name__ == "__main__":
    resolutions = [int(r) for r in sys.argv[1:]] or [200, 400, 800]
    for res in resolutions:
        run(res)
