# 3D FDTD (Meep) vs RCWA for the 8 unit cells

Phase: relative to the empty cell; the common offset between the methods (circular mean of the differences) is removed, because only the phase steps between levels affect the vortex plate.

## Resolution 40 px/µm (common offset -2.3°)

| Level | d (nm) | T0 FDTD | T0 RCWA | Phase difference after offset (°) |
|---|---|---|---|---|
| 0 | 310 | 0.997 | 0.994 | +5.5 |
| 1 | 140 | 0.990 | 0.993 | -2.6 |
| 2 | 168 | 0.999 | 0.996 | -4.8 |
| 3 | 188 | 0.981 | 0.976 | -3.4 |
| 4 | 200 | 0.974 | 0.973 | -3.2 |
| 5 | 226 | 0.983 | 0.984 | +3.0 |
| 6 | 244 | 0.969 | 0.958 | +1.3 |
| 7 | 284 | 0.991 | 0.979 | +4.2 |

RMS phase difference 3.7°, max |ΔT0| 0.012

## Resolution 80 px/µm (common offset -5.7°)

| Level | d (nm) | T0 FDTD | T0 RCWA | Phase difference after offset (°) |
|---|---|---|---|---|
| 0 | 310 | 0.995 | 0.994 | +1.7 |
| 1 | 140 | 0.990 | 0.993 | +0.2 |
| 2 | 168 | 0.999 | 0.996 | -3.0 |
| 3 | 188 | 0.979 | 0.976 | -2.1 |
| 4 | 200 | 0.972 | 0.973 | -1.8 |
| 5 | 226 | 0.986 | 0.984 | +1.4 |
| 6 | 244 | 0.967 | 0.958 | +1.9 |
| 7 | 284 | 0.988 | 0.979 | +1.6 |

RMS phase difference 1.9°, max |ΔT0| 0.009
