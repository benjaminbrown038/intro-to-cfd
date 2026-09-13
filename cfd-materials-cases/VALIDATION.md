# Validation record

All five default cases and all 32 comparison/refinement runs completed successfully.
13 physics tests passed. See `verification.log` for the individual tests.
Tested environment: Python, NumPy, SciPy, and Matplotlib versions are saved in each
result JSON. No claim of verification on other environments is made.

## Default numerical checks

| Check | Measured value |
|---|---:|
| Channel flow relative error | 0.00125 |
| Duct linear residual | 1.59531e-14 |
| Porous mass relative imbalance | 6.82811e-13 |
| Thermal energy relative imbalance | 1.232e-13 |
| Compliant tube maximum hoop strain | 0.025 |

Channel flow error decreases by a factor of four when cell count doubles, consistent
with second-order convergence. Square-duct flow is checked against the analytical
series, with <0.2% error at 64×64. Porous uniform and layered media match independent
linear-pressure and series-resistance references. The thermal tests check global
energy conservation, temperature bounds, uniform-temperature equilibrium, and a
conductivity trend. The compliant case matches an independent exact polynomial
integral and approaches the rigid-tube limit.

## Mesh limitations

The masked insert study changes flow by about 4.2% between 32×32 and 64×64; the porous
inclusion study changes flow by about 2.3% between 40×40 and 80×80. These curved-mask
cases are **not established as mesh independent**. Staircase geometry errors can be
nonmonotonic. Refine further or introduce boundary-fitted/cut-cell geometry before
using precise performance comparisons.

The thermal refinement study approaches outlet temperature from below (about
317.139, 317.241, and 317.299 K). This is evidence of refinement sensitivity, not an
independent validation of a real heated channel. Upwind numerical diffusion remains.

The compliant solver is checked against its own reduced governing model through an
independent analytic integral; this does not validate that model against experiments
or a full structural/fluid solver. No case has experimental validation, turbulence,
CAD meshing, transient flow, or a general moving-wall solver.
