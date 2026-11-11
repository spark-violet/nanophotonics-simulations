"""Export the metalens layout to GDS with gdsfactory.

Reads results/pillars.csv (written by design_metalens.py) and places one
circle per pillar. Layers:
  (1, 0)  pillars (TiO2), to be written by e-beam lithography
  (2, 0)  lens outline (for inspection only)
  (3, 0)  label text and alignment crosses

Outputs: results/metalens.gds and results/metalens_layout.png
"""
from pathlib import Path

import gdsfactory as gf
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle

import design

RES = Path(__file__).resolve().parents[1] / "results"
LAYER_PILLAR = (1, 0)
LAYER_OUTLINE = (2, 0)
LAYER_MARKS = (3, 0)
CIRCLE_POINTS = 32      # vertices per pillar polygon


def build():
    data = np.loadtxt(RES / "pillars.csv", delimiter=",", skiprows=1)
    gf.gpdk.PDK.activate()
    c = gf.Component("metalens_TiO2_f40um_D40um")

    # one cell per distinct diameter, then references -> small GDS file
    cells = {}
    for d_nm in np.unique(np.round(data[:, 2], 1)):
        cell = gf.Component(f"pillar_d{d_nm:.1f}nm".replace(".", "p"))
        cell.add_polygon(_circle(d_nm / 2e3), layer=LAYER_PILLAR)
        cells[d_nm] = cell
    for x, y, d_nm in data:
        ref = c.add_ref(cells[round(d_nm, 1)])
        ref.move((x, y))

    R = design.LENS_DIAMETER / 2
    c.add_polygon(_ring(R + 0.5, R + 1.0), layer=LAYER_OUTLINE)

    # alignment crosses and label outside the lens
    for sx, sy in [(-1, -1), (1, -1), (-1, 1), (1, 1)]:
        cross = gf.components.cross(length=4, width=0.5, layer=LAYER_MARKS)
        c.add_ref(cross).move((sx * (R + 6), sy * (R + 6)))
    label = gf.components.text(
        f"f={design.FOCAL_LENGTH:.0f}um D={design.LENS_DIAMETER:.0f}um "
        f"lambda={design.WAVELENGTH*1e3:.0f}nm", size=2, layer=LAYER_MARKS)
    c.add_ref(label).move((-R, -R - 10))

    path = RES / "metalens.gds"
    c.write_gds(path)
    return data, path


def _circle(radius):
    t = np.linspace(0, 2 * np.pi, CIRCLE_POINTS, endpoint=False)
    return list(zip(radius * np.cos(t), radius * np.sin(t)))


def _ring(r_in, r_out, n=256):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    outer = list(zip(r_out * np.cos(t), r_out * np.sin(t)))
    inner = list(zip(r_in * np.cos(t[::-1]), r_in * np.sin(t[::-1])))
    return outer + [outer[0]] + inner + [inner[0]]


def preview(data):
    """Full lens and a zoom on the centre, drawn from the pillar list."""
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.4), dpi=150)
    for ax, half in zip(axes, [design.LENS_DIAMETER / 2 + 1, 2.0]):
        sel = (np.abs(data[:, 0]) <= half + 0.5) & (np.abs(data[:, 1]) <= half + 0.5)
        patches = [Circle((x, y), d / 2e3) for x, y, d in data[sel]]
        ax.add_collection(PatchCollection(patches, facecolor="#7a1530", edgecolor="none"))
        ax.set_xlim(-half, half); ax.set_ylim(-half, half)
        ax.set_aspect("equal")
        ax.set_xlabel("x (µm)"); ax.set_ylabel("y (µm)")
    axes[0].set_title(f"Layout: {len(data)} pillars")
    axes[1].set_title("Centre (4 µm × 4 µm)")
    fig.tight_layout()
    fig.savefig(RES / "metalens_layout.png")


if __name__ == "__main__":
    data, path = build()
    preview(data)
    print(f"wrote {path.name} ({path.stat().st_size/1e3:.0f} kB), {len(data)} pillars")
