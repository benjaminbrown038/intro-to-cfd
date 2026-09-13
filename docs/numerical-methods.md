# Numerical methods and scope

## Transport (Tiers 1–2)

Tier 1 uses forward Euler with centered second differences for u_t=alpha*u_xx. Tier 2 uses conservative face-flux differences for u_t+a*u_x=alpha*u_xx, with periodic boundaries. Upwind and Lax–Wendroff advective fluxes are compared. Both add an explicit diffusive flux. See each tier's analytical reference.

## Flow (Tiers 3–5)

Equations: -mu*Laplacian(u) + grad(p) = f, and div(u)=0.

This is steady incompressible creeping flow. Velocity has two components in 2D and three in 3D. Pressure is cell-centered; each velocity component occupies the internal faces normal to its own coordinate axis. Pressure is fixed to zero at the first cell to remove its arbitrary constant.

The integrated divergence matrix B uses a shared oriented face flux for the two neighboring cells. With velocity diffusion block K and prescribed boundary flux contribution q, the saddle system is:

K*u - B^T*p = f; B*u = q.

One redundant continuity equation and its pressure unknown are removed. Velocity diffusion conductances are mu*area/distance on dual control volumes. Tangential wall velocities contribute through the physical wall distance. The result reports divergence in inverse seconds and the normalized algebraic residual. Small residuals do not by themselves establish small discretization error.

SciPy sparse direct solves are used in Python; MATLAB uses sparse backslash. Python reference: [SciPy spsolve](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.spsolve.html).

## Adaptation

Velocity-gradient activity is summed over components and directions. For each axis, volume-integrated activity in an interval is multiplied by that interval's width squared. The highest-scoring intervals receive midpoint planes. Grids stay tensor-product and conforming. Solving each new steady problem from scratch avoids any requirement for transferring transient state. This is not octree or block AMR and the activity score is not a rigorous error estimator.

## Obstacle

Tier 4 adds penalty*mask*u to momentum inside a circular region. This is a porous-resistance approximation. Study both spatial resolution and penalty strength before interpreting surface forces. The example does not model an unsteady wake.
