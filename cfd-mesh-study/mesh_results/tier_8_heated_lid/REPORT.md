# Mesh study: tier_8_heated_lid

Refine further: one or more finest-pair changes exceed the threshold or are undefined.

| Mesh | Final status | Kinetic energy | Hot-wall Nu | Velocity change from preceding mesh |
|---|---|---:|---:|---:|
| 16 × 16 | steady_tolerance_reached | 0.038287 | 2.945256 | — |
| 32 × 32 | steady_tolerance_reached | 0.043032 | 2.919365 | 7.005% |
| 48 × 48 | steady_tolerance_reached | 0.044581 | 2.915517 | 2.460% |

Common dimensionless time-step cap: 0.00112484. Change threshold: 1%.

Pairwise field differences use fine-grid interpolation onto coarse centers. They are not errors against an exact solution. Runtime depends on grid cost and the number of steps needed to reach the stopping criterion.
