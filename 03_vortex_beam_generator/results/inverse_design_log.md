# Inverse design log

Gradient check (one pixel): autograd 8.331083e-03, finite difference 8.331084e-03, relative difference 1.2e-07

Score = sqrt(T0)·cos(phase error), 1 is perfect. 'Worst ±5 nm' = lowest score when every edge moves 5 nm in or out (pillar diameter ±10 nm).

| Level | Pillar score | Pillar worst ±5 nm | Freeform score | Freeform T0 | Freeform phase error (°) | Freeform worst ±5 nm |
|---|---|---|---|---|---|---|
| 0 | 0.997 | 0.924 | 0.992 | 0.992 | +4.9 | 0.867 |
| 1 | 0.996 | 0.972 | 0.947 | 0.983 | -17.2 | 0.889 |
| 2 | 0.998 | 0.915 | 0.979 | 0.999 | -11.5 | 0.863 |
| 3 | 0.988 | 0.891 | 0.972 | 0.980 | -11.0 | 0.803 |
| 4 | 0.966 | 0.740 | 0.984 | 0.972 | +3.5 | 0.881 |
| 5 | 0.988 | 0.872 | 0.993 | 0.987 | -1.8 | 0.829 |
| 6 | 0.978 | 0.901 | 0.964 | 0.949 | +8.3 | 0.867 |
| 7 | 0.977 | 0.840 | 0.991 | 0.989 | -4.7 | 0.945 |

8-level FOM (zeroth-order efficiency of the phase plate): pillars 0.972, freeform 0.960
