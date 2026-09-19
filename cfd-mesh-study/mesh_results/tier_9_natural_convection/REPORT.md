# Mesh study: tier_9_natural_convection

Refine further: one or more finest-pair changes exceed the threshold or are undefined.

| Mesh | Final status | Kinetic energy | Hot-wall Nu | Velocity change from preceding mesh |
|---|---|---:|---:|---:|
| 16 × 16 | steady_tolerance_reached | 0.008816 | 2.328465 | — |
| 32 × 32 | steady_tolerance_reached | 0.008913 | 2.265555 | 4.510% |
| 48 × 48 | steady_tolerance_reached | 0.008977 | 2.253709 | 1.158% |

Common dimensionless time-step cap: 0.00127086. Change threshold: 1%.

Pairwise field differences use fine-grid interpolation onto coarse centers. They are not errors against an exact solution. Runtime depends on grid cost and the number of steps needed to reach the stopping criterion.
