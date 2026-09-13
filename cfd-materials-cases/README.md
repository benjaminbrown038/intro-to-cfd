# Intro to CFD — Materials and Geometry Cases

Five runnable Python cases for `intro-to-cfd`, progressing from a verified channel
solver to heterogeneous porous flow, conjugate heat transfer, and reduced
fluid–wall coupling. These are steady educational solvers with actual numerical
solutions, editable JSON inputs, plots, data exports, and physics tests.

## Add to your repository

Extract `cfd-materials-cases.zip`, then place the complete `cfd-materials-cases`
folder inside your local `intro-to-cfd` repository. This folder is self-contained;
it does not require replacing your existing README, requirements, or tier folders.
The existing repository files were not available when this pack was built, so no
claim is made that it has been merged or pushed. This addition is Python only.

From a terminal inside the extracted folder:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m cfd_cases --all
python -m unittest discover -s tests -v
```

On Windows, activate using `.venv\Scripts\activate` instead. A standard CPU is
sufficient for the supplied cases. No GPU, OpenFOAM, or commercial solver is needed.

Run an individual case or all comparison studies:

```bash
python -m cfd_cases --case cases/tier-4/case.json
python -m cfd_cases --studies
python -m cfd_cases --all --out my_results
```

Commands must run from this folder. Results with the same name in the same output
folder are replaced on rerun; use `--out` to preserve a separate experiment.
To save a new experiment, copy a case JSON, change its `name` and inputs, and pass
its path to `--case`. Do not change the `model` unless switching solvers.

## Contents and progression

| Tier | Configuration | Model | Material / geometry choices |
|---|---|---|---|
| 1 | `cases/tier-1/case.json` | 1D finite-volume parallel-plate flow | Viscosity, height, pressure gradient |
| 2 | `cases/tier-2/case.json` | 2D cross-sectional duct velocity | Rectangle, ellipse, longitudinal cylindrical insert |
| 3 | `cases/tier-3/case.json` | 2D heterogeneous Darcy flow | Uniform medium, low-permeability layer, elliptical inclusion |
| 4 | `cases/tier-4/case.json` | 2D solid–fluid energy equation | Solid conductivity, wall thickness, fluid thermal properties |
| 5 | `cases/tier-5/case.json` | Reduced 1D compliant tube | Elastic modulus, thickness, radius, pressure |

- `cfd_cases/fv.py`: shared conservative diffusion assembly and face fluxes.
- `cfd_cases/models.py`: governing models and physical input checks.
- `cfd_cases/__main__.py`: command line, plotting, exports, parameter studies.
- `tests/test_physics.py`: analytical, conservation, limiting-case, and refinement checks.
- `results/`: precomputed examples and studies, included for immediate inspection;
  ignored by git by default because they can be regenerated.
- `VALIDATION.md`: measured checks and remaining numerical limitations.

## Tier 1 — Pressure-driven channel

Solve `-mu d²u/dy² = G`, where `G = -dp/dx > 0`, with zero velocity on both
stationary plates. The channel is infinite in the out-of-plane direction and
fully developed along the flow direction. The exact reference is
`u(y) = G y (H-y)/(2 mu)` and flow per unit depth is `G H³/(12 mu)`.

The cell-centered solver balances viscous fluxes against the pressure driving
term. Boundary faces are half a cell from their adjacent centers. Summed
velocity times cell height gives flow **per unit depth**, in m²/s.

Try doubling viscosity, halving height, or doubling pressure gradient. Keep the
flow laminar; the reported hydraulic Reynolds number uses `Dh = 2H`. Density
changes this diagnostic but does not change this fixed-gradient viscous solution.

## Tier 2 — Cross-section geometry

Solve `-mu (d²u/dx² + d²u/dy²) = G` across the duct cross-section. Here `u` is
velocity **along the duct**, perpendicular to the plotted plane. No-slip walls
bound all fluid cells. An insert is a solid cylinder extending along the whole
duct, not an obstacle facing an inlet. This model does not compute wakes, bends,
streamwise contractions, entrance effects, or recirculation.

Set `shape` to `rectangle`, `ellipse`, or `insert`. Edit `width`, `height`, and
`insert_radius`. All geometry comparisons use the same bounding dimensions and
pressure gradient; their open areas differ. Conductance differences therefore
include both area and shape effects. Flow is in m³/s.

Curved boundaries use a Cartesian staircase mask. Refinement changes both the
numerical solution and represented shape. The square-duct test uses an independent
Fourier-series solution. Curved cases require separate geometric convergence.
Changing a solid material name alone would have no effect in this rigid,
unheated, no-slip model; viscosity and geometry are the active inputs.

## Tier 3 — Permeability and flow diversion

Solve `div(v) = 0`, `v = -(K/mu) grad(p)` in a fully saturated porous rectangle.
Left pressure is prescribed, right pressure is zero gauge, and top/bottom have
zero normal flux. All cells are porous material; this is not a free-fluid flow
solver with a filter inserted into an otherwise empty channel.

`permeability` sets the background K in m². `contrast` multiplies K inside an
inclusion or barrier. For example, 0.01 gives a material 100 times less permeable.
The barrier spans x/L = 0.4–0.6; the inclusion is centered with semi-axes 0.15L and
0.30H. Face mobility uses the harmonic mean, representing resistances in series.

Saved velocity is **Darcy velocity / superficial velocity**, not pore velocity.
Porosity is not an input because this steady model solves superficial discharge;
pore velocity would also require porosity. No gravity, inertia, anisotropy, or
Forchheimer correction is modeled. Face-based conservation metrics are calculated
before interpolating velocities to cell centers for streamlines.

## Tier 4 — Solid–fluid heat transfer

Solve `rho cp u dT/dx = div(k grad(T))` in the fluid and conduction in two solid
slabs. Conductivity changes by material. A single finite-volume system enforces
perfect thermal contact; harmonic face conductivity keeps interface heat flux
consistent. Upwind advection is first order. Axial and transverse conduction are
included. The fluid velocity is a prescribed, discretely normalized, fully
developed parabolic profile.

Boundary conditions:

- Fluid inlet: prescribed temperature, including advective and conductive flux.
- Fluid outlet: zero conductive flux; advective energy leaves the domain.
- Solid inlet/outlet ends: insulated.
- Outer top and bottom surfaces: prescribed temperature.

The energy check includes **outer-wall heat + inlet conduction = advective gain**.
Ignoring inlet conduction would give an incorrect apparent imbalance.

Compare `k_solid = 15, 205, 400 W/(m K)`, illustrative values broadly representing
stainless-steel-like, aluminum-like, and copper-like conductivities. These are
idealized constant inputs, not a grade- or temperature-specific material database.
The high-conductivity cases may differ only slightly when fluid-side resistance
dominates. Material properties stay constant: temperature does not feed back into
viscosity or momentum in this version.

The grid must align with both interfaces:
`fluid_height / n_fluid = wall_thickness / n_wall`.
To double wall thickness without changing cell spacing, double `n_wall` too.
For refinement multiply `nx`, `n_fluid`, and `n_wall` by the same factor.

## Tier 5 — Pressure and compliant walls

A thin membrane tube expands according to
`R(p) = R0 [1 + (p-p_external) R0/(E t)]`.
Local laminar resistance satisfies `dp/dx = -8 mu Q/(pi R(p)^4)`.
Integrating along the tube gives
`Q = pi/(8 mu L) integral[p_out,p_in] R(p)^4 dp`.
The code performs trapezoidal quadrature, reconstructs position from cumulative
resistance, and interpolates pressure and radius onto an axial output grid.

Both pressures are prescribed, so expansion increases Q relative to a rigid tube.
The model assumes zero axial membrane stress, small hoop strain, negligible
inertia, slowly varying geometry, constant thickness to first order, and local
Poiseuille behavior. It has no wall dynamics, wave propagation, collapse,
viscoelasticity, axial wall coupling, or moving mesh. It is a **reduced steady
coupling example**, not a general fluid–structure interaction solver.

Input guards reject compressive transmural pressure, thickness/radius > 0.1,
hoop strain > 5%, and maximum Reynolds number > 1000. These guards do not establish
validity for every user-defined case; retain a long slender tube and gentle radius
variation. Use `young_modulus` to compare stiffness. A very large modulus recovers
the rigid-tube result. Default fluid viscosity is 0.01 Pa s, not water viscosity.

## Outputs and experiments

Every run writes:

| File | Contents |
|---|---|
| `result.json` | Complete case settings, metrics, runtime, environment versions |
| `fields.csv` | Named SI-unit columns, flattened row-major for 2D fields |
| `fields.npz` | Arrays with original shapes for later analysis |
| `plot.png` | Velocity, pressure, materials, temperature, or wall response |

`--studies` runs 32 cases across 10 studies: viscosity, cross-section shape,
permeability, wall conductivity, stiffness, and resolution for each solver. It also
writes `studies.csv` and `study_comparisons.png`. These are numerical experiments;
small algebraic residuals and conservation errors alone do not prove mesh accuracy.

```python
import numpy as np
fields = np.load('results/tier_4_heated_wall/fields.npz')
temperature = fields['temperature_K']
print(temperature.shape, temperature.min(), temperature.max())
```

## Numerical references and next extensions

The finite-volume formulation follows the control-volume balance of face diffusion
and convection described in [NIST FiPy's finite-volume documentation](https://pages.nist.gov/fipy/en/stable/numerical/discret.html).
This implementation uses NumPy/SciPy directly; FiPy is not a dependency.
The porous model uses the [Darcy pressure–velocity formulation](https://doc.comsol.com/6.3/doc/com.comsol.help.porous/porous_ug_fluidflow_porous.07.004.html).

A future bend/obstacle case needs a pressure–velocity Navier–Stokes solver; surface
roughness in turbulent pipes needs a turbulence/wall model; transient flexible
walls need structural dynamics and fluid–structure coupling. These are separate
extensions, not hidden capabilities of the current tiers. This pack provides
verified building blocks and explicitly reduced multiphysics cases to develop
before those steps.
