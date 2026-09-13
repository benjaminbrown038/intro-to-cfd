# Increasing complexity

Every tier now has a runnable example. Each tier's broader progression remains incremental.

| Tier | Implemented | Next extension |
| --- | --- | --- |
| 1 | Transient diffusion with analytical checks | Steady sources and independent time refinement |
| 2 | Periodic transport; two advective fluxes | Limiters and discontinuous-front cases |
| 3 | Coupled 2D steady Stokes cavity | Momentum advection and time integration |
| 4 | Analytical channel check; penalized circular obstacle | Fitted geometry, developing flow, unsteady wakes |
| 5 | 2D/3D adaptive tensor meshes; uniform comparisons | Local block/octree refinement, conservative transient transfer |

Increase physical complexity, geometry complexity, numerical complexity, and computational cost separately. Turbulence, compressibility, multiphase flow, fluid–structure interaction, GPU acceleration, and distributed computing are future extensions, not implemented features.
