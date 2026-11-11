# Vortex-beam generator: results

Plate: 30 µm diameter, focal length 30 µm, 8 phase levels, TiO2 cells H = 650 nm, P = 400 nm, λ = 633 nm.

| Unit cells | l | OAM purity | Efficiency (r < 3 µm) |
|---|---|---|---|
| ideal_continuous | 1 | 1.000 | 0.826 |
| ideal_8_level | 1 | 0.999 | 0.784 |
| pillars_rcwa | 1 | 0.999 | 0.762 |
| pillars_fdtd | 1 | 0.999 | 0.761 |
| freeform_inverse | 1 | 0.998 | 0.748 |
| ideal_continuous | 2 | 1.000 | 0.818 |
| ideal_8_level | 2 | 1.000 | 0.777 |
| pillars_rcwa | 2 | 1.000 | 0.755 |
| pillars_fdtd | 2 | 0.999 | 0.755 |
| freeform_inverse | 2 | 0.999 | 0.741 |

## Fabrication tolerance (l = 1)

Every edge moved by the same amount (systematic over- or under-etch).

| Unit cells | Edge bias | OAM purity | Efficiency (r < 3 µm) |
|---|---|---|---|
| pillars_rcwa | -5 nm | 0.999 | 0.747 |
| freeform_inverse | -5 nm | 0.997 | 0.694 |
| pillars_rcwa | +5 nm | 0.999 | 0.758 |
| freeform_inverse | +5 nm | 0.998 | 0.735 |

OAM purity: share of the focal-plane power in the target mode (azimuthal decomposition on rings, r < 3 µm). Efficiency: power within r < 3 µm divided by the power incident on the plate.
