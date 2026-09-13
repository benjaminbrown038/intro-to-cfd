# Verification and measurements

A run finishing without an error is not enough to mark a solver complete.

- Check units and prescribed boundary conditions.
- Compare against an analytical solution where available.
- Measure mass conservation and velocity divergence for incompressible flow.
- Refine the mesh and time step separately to assess numerical error.
- Distinguish equation residuals from error relative to a reference solution.
- Record assembly, pressure-solve, time-step, and total run timings separately when relevant.
- For adaptive meshes, verify conservation across refinement interfaces and compare accuracy against uniform refinement.
- Keep physical assumptions and Reynolds number explicit; convergence alone does not validate those assumptions.

These are educational implementations, not validated engineering design tools.
