# Increasing complexity

Build and verify each step before moving to the next. All steps below are planned.

| Tier | Progression within the tier | Completion evidence |
| --- | --- | --- |
| 1 | Steady 1D diffusion → transient diffusion → mesh/time-step studies | Analytical agreement and measured error trends |
| 2 | Pure advection → advection–diffusion → scheme comparisons | Transport accuracy, conservation, and stability measurements |
| 3 | Pressure equation exercise → low-Re cavity flow → grid refinement | Small divergence, converged residuals, and a documented benchmark comparison |
| 4 | Fully developed channel → developing channel → cylinder flow → unsteady wake where physically appropriate | Channel analytical agreement, mass balance, and converged force/time histories |
| 5 | 2D uniform-grid baseline → local refinement → conservative adaptive cycle → 3D uniform baseline → 3D adaptation | Interface conservation and accuracy/cost comparison against uniform refinement |

## Increase one difficulty at a time

Physical complexity: diffusion, transport, incompressible flow, unsteady flow.
Geometric complexity: line, rectangle, obstacle, three-dimensional domain.
Numerical complexity: uniform grids, coupled equations, local refinement, adaptive cycles.
Computational complexity: measure cost first; optimize only observed bottlenecks.

Turbulence models, compressible flow, multiphase flow, fluid–structure interaction, GPU acceleration, and distributed computing are future extensions. They are not prerequisites for the five tiers and do not have placeholder directories in this starter.
