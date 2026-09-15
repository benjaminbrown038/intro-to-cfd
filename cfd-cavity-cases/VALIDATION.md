# Numerical validation record

Three default cases and five verification-study runs completed. Nine automated
unit/physics tests passed. This is educational numerical verification, not
validation of an aircraft component or a general CFD solver.

## Natural-convection benchmark

Pr=0.71, Ra=1000, stationary no-slip walls, hot left/cold right, insulated top/bottom.
Reference mean wall Nusselt number: 1.118, as reported in the
[FEATool developer tutorial](https://www.featool.com/doc/Multiphysics_04_natural_convection1)
for the [De Vahl Davis benchmark](https://doi.org/10.1002/fld.1650030305).

| Grid | Computed Nu | Relative error | Stopping status |
|---|---:|---:|---|
| 16 × 16 | 1.148519 | 2.730% | steady_tolerance_reached |
| 32 × 32 | 1.130902 | 1.154% | steady_tolerance_reached |
| 48 × 48 | 1.126036 | 0.719% | steady_tolerance_reached |

Error decreases under refinement. The 48×48 grid is within 1% of the reported
reference heat-transfer value. This single integral quantity does not establish
accuracy for all velocities, pressures, parameters, or transient behavior.
The Ra=10000 default case gives Nu≈2.2656, about 1.01% above the reported 2.243.

## Default termination status

| Case | Status | Final field-change rate |
|---|---|---:|
| tier_6_lid_cavity | steady_tolerance_reached | 9.43e-07 |
| tier_8_heated_lid | end_time_reached | 3.3e-06 |
| tier_9_natural_convection | steady_tolerance_reached | 9.11e-07 |

The heated moving-lid run completes t=40 but has not met the 1e-6 steady criterion.
It is saved as a completed transient, not a converged steady solution. Increase
end_time to continue studying its approach to steady flow; reruns start from rest.

## Checks and limits

- Projection removes a manufactured gradient field and reduces random-field
  divergence below 1e-9 while keeping normal wall velocities zero.
- Quiescent equilibrium is preserved.
- Linear conduction is preserved, with Nu=1.
- A sinusoidal conduction eigenmode shows second-order spatial operator convergence.
- Conservative thermal fluxes close the per-step stored-energy balance below 1e-12.
- The hot side rises and the cold side descends in a buoyancy startup test.
- Invalid inputs and incomplete max-step termination are checked.
- Default maximum discrete divergence is below 1e-13.

At Re=100, n=24, t=5, halving cfl from 0.4 to 0.2 changes integrated kinetic energy
by about 0.005%. This is one temporal sensitivity check, not proof of temporal
accuracy in other regimes. Time is nondimensional and this comparison is transient.

Momentum advection is first order and dissipative; only diffusion has a demonstrated
second-order spatial operator test. No complete Ghia centerline benchmark has been
performed. Pressure, stresses, and lid-corner behavior are not independently
validated. The Boussinesq approximation assumes small density variations and constant
transport properties. High Mach number, turbulence, FSI, solid conduction, variable
viscosity, and boiling are not represented.

Per-run result JSON files record the actual software environment and runtimes.
