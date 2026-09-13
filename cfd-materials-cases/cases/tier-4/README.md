# Tier 4: Heated solid walls

Change k_solid to 15, 205, or 400. Inspect the temperature field and energy balance. Keep wall_thickness/n_wall equal to fluid_height/n_fluid. Velocity is prescribed; temperature does not change viscosity.

From the package root:

```bash
python -m cfd_cases --case cases/tier-4/case.json
```

See the root README for governing equations, units, boundary conditions, and assumptions.
