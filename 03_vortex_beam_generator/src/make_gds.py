"""GDS layouts of the vortex-beam generator (l = 1).

Two files, one per unit-cell set:
  vortex_pillars.gds   8 circular pillar cells
  vortex_freeform.gds  8 inverse-designed cells; their pixel masks are
                       merged into clean polygons with KLayout before export

Each distinct unit cell is stored once and placed by reference, which keeps
the files small. Layers: (1,0) TiO2, (2,0) plate outline, (3,0) marks/label.

Outputs (in ../results): the two GDS files and vortex_layout.png
"""
from pathlib import Path

import gdsfactory as gf
import klayout.db as kdb
import matplotlib.pyplot as plt
import numpy as np

import config as C
from pillar_library import load_geometry
from vortex import plate_phase

RES = Path(__file__).resolve().parents[1] / "results"
L_TIO2, L_OUTLINE, L_MARK = (1, 0), (2, 0), (3, 0)


def lattice_levels(P, charge=C.TOPOLOGICAL_CHARGE):
    n = int(np.floor(C.PLATE_DIAMETER / 2 / P))
    idx = np.arange(-n, n + 1) * P
    X, Y = np.meshgrid(idx, idx, indexing="ij")
    inside = X**2 + Y**2 <= (C.PLATE_DIAMETER / 2) ** 2
    phi = plate_phase(X, Y, charge)
    level = np.mod(np.round(phi / (2 * np.pi / C.N_LEVELS)), C.N_LEVELS).astype(int)
    return X[inside], Y[inside], level[inside]


def circle(r, n=32):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return list(zip(r * np.cos(t), r * np.sin(t)))


def mask_polygons(mask, pixel, P):
    """Merge a binary pixel mask into polygons (cell centred at the origin)."""
    dbu = 0.001
    region = kdb.Region()
    n = mask.shape[0]
    for i, j in zip(*np.nonzero(mask)):
        x0 = -P / 2 + i * pixel
        y0 = -P / 2 + j * pixel
        region.insert(kdb.Box(round(x0 / dbu), round(y0 / dbu),
                              round((x0 + pixel) / dbu), round((y0 + pixel) / dbu)))
    region.merge()
    polys = []
    for poly in region.each():
        hull = [(p.x * dbu, p.y * dbu) for p in poly.each_point_hull()]
        holes = [[(p.x * dbu, p.y * dbu) for p in poly.each_point_hole(h)]
                 for h in range(poly.holes())]
        polys.append((hull, holes))
    return polys


def frame(c, name):
    R = C.PLATE_DIAMETER / 2
    t = np.linspace(0, 2 * np.pi, 256, endpoint=False)
    outer = list(zip((R + 1) * np.cos(t), (R + 1) * np.sin(t)))
    inner = list(zip((R + 0.5) * np.cos(t[::-1]), (R + 0.5) * np.sin(t[::-1])))
    c.add_polygon(outer + [outer[0]] + inner + [inner[0]], layer=L_OUTLINE)
    for sx, sy in [(-1, -1), (1, -1), (-1, 1), (1, 1)]:
        c.add_ref(gf.components.cross(length=4, width=0.5, layer=L_MARK)).move(
            (sx * (R + 6), sy * (R + 6)))
    c.add_ref(gf.components.text(name, size=2, layer=L_MARK)).move((-R, -R - 10))


def build(kind, H, P, X, Y, level):
    c = gf.Component(f"vortex_l{C.TOPOLOGICAL_CHARGE}_{kind}")
    cells = []
    if kind == "pillars":
        d = np.loadtxt(RES / "levels_pillars.csv", delimiter=",", skiprows=1)[:, 2] / 1e3
        for k in range(C.N_LEVELS):
            cell = gf.Component(f"pillar_level{k}")
            cell.add_polygon(circle(d[k] / 2), layer=L_TIO2)
            cells.append(cell)
        shapes = [("circle", dk / 2) for dk in d]
    else:
        z = np.load(RES / "freeform_cells.npz")
        shapes = []
        for k, rho in enumerate(z["designs"]):
            cell = gf.Component(f"freeform_level{k}")
            polys = mask_polygons(rho.astype(bool), float(z["pixel"]), P)
            for hull, holes in polys:
                if holes:
                    # gdsfactory polygons have no holes: subtract via boolean
                    tmp = gf.Component()
                    tmp.add_polygon(hull, layer=L_TIO2)
                    for h in holes:
                        hc = gf.Component()
                        hc.add_polygon(h, layer=L_TIO2)
                        tmp = gf.boolean(tmp, hc, "not", layer=L_TIO2)
                    cell.add_ref(tmp)
                else:
                    cell.add_polygon(hull, layer=L_TIO2)
            cells.append(cell)
            shapes.append(("mask", rho.astype(bool)))
    for x, y, k in zip(X, Y, level):
        c.add_ref(cells[k]).move((x, y))
    frame(c, f"vortex l={C.TOPOLOGICAL_CHARGE} f={C.FOCAL_LENGTH:.0f}um {kind}")
    path = RES / f"vortex_{kind}.gds"
    c.write_gds(path)
    return path, shapes


def check_gds(path, expected_cells):
    """Read the file back and count TiO2 cell placements."""
    ly = kdb.Layout()
    ly.read(str(path))
    li = ly.find_layer(*L_TIO2)
    count = sum(1 for inst in ly.top_cell().each_inst()
                if not ly.cell(inst.cell_index).bbox_per_layer(li).empty())
    return count, expected_cells


def raster_cell(shape, P, npix=40):
    typ, val = shape
    if typ == "mask":
        return val.astype(float)
    x = (np.arange(npix) + 0.5) / npix * P - P / 2
    XX, YY = np.meshgrid(x, x, indexing="ij")
    return (XX**2 + YY**2 <= val**2).astype(float)


def preview(H, P, X, Y, level, shapes_by_kind):
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.2), dpi=150)
    R = C.PLATE_DIAMETER / 2
    n = int(np.floor(R / P))
    grid = np.full((2 * n + 1, 2 * n + 1), np.nan)
    grid[np.round(X / P).astype(int) + n, np.round(Y / P).astype(int) + n] = level
    im = axes[0].imshow(grid.T, origin="lower", cmap="twilight", interpolation="nearest",
                        extent=[-(n + .5) * P, (n + .5) * P, -(n + .5) * P, (n + .5) * P])
    axes[0].set_title(f"Phase level per cell, l = {C.TOPOLOGICAL_CHARGE}")
    fig.colorbar(im, ax=axes[0], fraction=0.046, label="level")
    m = int(round(2.0 / P))                       # cells on each side of the centre
    for ax, (kind, shapes) in zip(axes[1:], shapes_by_kind.items()):
        tiles = {k: raster_cell(shapes[k], P) for k in range(C.N_LEVELS)}
        npx = tiles[0].shape[0]
        img = np.zeros(((2 * m + 1) * npx, (2 * m + 1) * npx))
        for i in range(-m, m + 1):
            for j in range(-m, m + 1):
                g = grid[i + n, j + n]
                if not np.isnan(g):
                    img[(i + m) * npx:(i + m + 1) * npx, (j + m) * npx:(j + m + 1) * npx] = tiles[int(g)]
        half = (m + 0.5) * P
        ax.imshow(img.T, origin="lower", cmap="Reds", vmin=0, vmax=1.3,
                  extent=[-half, half, -half, half], interpolation="nearest")
        ax.set_title(f"{kind}: centre {2*half:.1f} µm × {2*half:.1f} µm")
    for ax in axes:
        ax.set_aspect("equal"); ax.set_xlabel("x (µm)"); ax.set_ylabel("y (µm)")
    fig.tight_layout()
    fig.savefig(RES / "vortex_layout.png")


def min_features(mask, pixel):
    """Smallest solid and gap width (nm) of a periodic binary cell.

    Uses morphological opening (solid) and closing (gap): the largest disk
    that leaves the pattern unchanged gives the narrowest feature."""
    from scipy.ndimage import binary_closing, binary_opening
    n = mask.shape[0]
    big = np.tile(mask, (3, 3))

    def widest(op, m):
        # largest disk for which the pattern is unchanged: the pattern is then a
        # union of such disks (solid) or its gaps are (void), so nothing in it
        # is narrower than the disk
        for r in range(n // 2, 0, -1):
            y, x = np.ogrid[-r:r + 1, -r:r + 1]
            disk = x**2 + y**2 <= r**2
            if np.array_equal(op(m, structure=disk)[n:2 * n, n:2 * n], m[n:2 * n, n:2 * n]):
                return (2 * r + 1) * pixel * 1e3
        return pixel * 1e3

    solid = widest(binary_opening, big) if mask.any() else np.nan
    gap = widest(binary_closing, big) if (~mask).any() else np.nan
    return solid, gap


def main():
    gf.gpdk.PDK.activate()
    H, P = load_geometry()
    X, Y, level = lattice_levels(P)
    shapes = {}
    for kind in ["pillars", "freeform"]:
        if kind == "freeform" and not (RES / "freeform_cells.npz").exists():
            continue
        path, shapes[kind] = build(kind, H, P, X, Y, level)
        n, exp = check_gds(path, len(X))
        print(f"{path.name}: {path.stat().st_size/1e3:.0f} kB, "
              f"{n} TiO2 cells placed (expected {exp})")
    preview(H, P, X, Y, level, shapes)

    if "freeform" in shapes:
        z = np.load(RES / "freeform_cells.npz")
        lines = ["# Manufacturability of the freeform cells", "",
                 "Narrowest solid feature and narrowest gap (periodic cell, 10 nm pixels; "
                 "a value of X nm means the shape is a union of X-nm disks). Fabrication limit: "
                 f"{C.MIN_FEATURE*1e3:.0f} nm features, {C.MIN_GAP*1e3:.0f} nm gaps.", "",
                 "| Level | Narrowest solid (nm) | Narrowest gap (nm) |", "|---|---|---|"]
        for k, rho in enumerate(z["designs"]):
            sol, gap = min_features(rho.astype(bool), float(z["pixel"]))
            lines.append(f"| {k} | {sol:.0f} | {gap:.0f} |")
        (RES / "manufacturability.md").write_text("\n".join(lines) + "\n")
        print("\n".join(lines))


if __name__ == "__main__":
    main()
