# Intro to CFD — Unsteady Cavity Cases

Three runnable 2D incompressible simulations that extend the materials/geometry
pack with pressure–velocity coupling, transient flow, heat advection, and buoyancy.
This is a separate self-contained folder inside the same `intro-to-cfd` repository.
It does not overwrite your existing `cfd-materials-cases` folder.

## Install and run

Extract the ZIP and put `cfd-cavity-cases` inside your local `intro-to-cfd` repo.
Open a terminal inside `cfd-cavity-cases`:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m cavity_cases --all --animate
```

On Windows, use `.venv\Scripts\activate` to activate the environment.
Python 3.9+ is supported by the source syntax. Dependencies must also support your
Python version; installed versions in the tested environment are recorded in JSON.
No GPU or commercial CFD package is required. The included examples use 32×32 cells.

```bash
python -m cavity_cases --case cases/tier-9/case.json --animate
python -m cavity_cases --studies --out verification_studies
python -m unittest discover -s tests -v
```

All paths above are relative to this folder. Copy a JSON and change its `name` to
make a new experiment. Use `--case path/to/new.json` to run it. Results with the same
name and output directory are overwritten on rerun; use `--out` for separate runs.

## Cases

| Tier | Case | Physical mechanism | Aerospace learning connection |
|---|---|---|---|
| 6 | Moving-lid cavity, Re=100 | Wall shear drives circulation | Foundation for internal viscous-flow solvers |
| 8 | Heated moving-lid cavity, Re=100, Pr=0.71, Ri=1 | Wall motion and buoyancy drive flow; flow transports heat | Interaction of forced circulation and heating |
| 9 | Natural convection, Ra=10000, Pr=0.71 | A hot wall drives rising fluid and circulation | Simplified enclosed equipment/air-cavity cooling |

Tier 7 remains reserved for the proposed step/cylinder solver with inlet/outlet
boundaries and more complex geometry. Tier 10's flexible wall is also not included.
The numbering follows our earlier roadmap. These are three implemented cases,
not a claim that every intervening tier is complete.

All geometries here are closed, rigid, unit-square cavities. The hot wall is left,
the cold wall is right, and the top/bottom are thermally insulated. The moving lid
travels to the right; the other walls are stationary. All walls are impermeable.
Tier 6 is isothermal. Thermal cases start with a linear hot-to-cold temperature
profile and zero velocity; they simulate flow startup, not sudden heating of an
initially uniform-temperature enclosure.

## What is new compared with the first pack?

The first duct case assumed fully developed axial flow. Its heated channel used
a prescribed velocity profile. Here the program advances both velocity components
and solves a pressure Poisson equation at every time step. In thermal cases,
velocity advects heat and temperature drives a Boussinesq buoyancy force.

No physical fluid or wall properties are temperature-dependent in this release.
Buoyancy approximates the effect of small density differences; density is otherwise
constant. This is not compressible gas flow, turbulence, variable-viscosity flow,
solid conduction, structural deformation, boiling, or general FSI. Hot bleed-air
and cryogenic applications need different/extended models.

## Equations and scales

All solver variables and JSON controls are dimensionless. These are not meters,
seconds, pascals, or kelvin until reference scales are supplied.

The solver advances

- `div(u) = 0`
- `du/dt + (u dot grad)u = -grad(p) + nu_star laplacian(u) + B (theta - 0.5) e_y`
- `dtheta/dt + div(u theta) = alpha_star laplacian(theta)`

where `theta = (T-T_cold)/(T_hot-T_cold)` and `p` is dimensionless kinematic
pressure with its spatial mean removed. The reference hydrostatic contribution
is absorbed into pressure. Pressure is a projection-method output with splitting
and boundary errors; it is not validated here as an accurate wall-load prediction.

| Mode | Reference speed | nu_star | alpha_star | B |
|---|---|---|---|---|
| lid | lid speed U | 1/Re | unused | 0 |
| thermal_lid | lid speed U | 1/Re | 1/(Re Pr) | Ri |
| natural | sqrt(g beta deltaT L) | sqrt(Pr/Ra) | 1/sqrt(Ra Pr) | 1 |

`Re = U L/nu`, `Pr = nu/alpha`, `Ra = g beta deltaT L³/(nu alpha)`, and
`Ri = g beta deltaT L/U²`. Here dimensional `nu = mu/rho` and
`alpha = k/(rho cp)` are diffusivities in m²/s; `beta` is the volumetric thermal
expansion coefficient in 1/K. Positive buoyancy is upward and gravity is downward.

To map outputs back to dimensions:

- length = dimensionless length × L
- velocity = dimensionless velocity × reference speed
- time = dimensionless time × L/reference speed
- temperature = T_cold + theta × deltaT
- dynamic pressure variation = dimensionless pressure × rho × reference speed²

Different fluids enter through Re, Pr, Ra, and Ri computed from actual properties.
For constant-property models, use properties appropriate to the operating range.
Boussinesq buoyancy requires small relative density changes, typically assessed
through `beta*deltaT << 1`, and incompressible modeling requires the relevant flow
conditions to justify it. There is no automatic physical applicability check because
raw dimensional properties are not inputs to this solver.

## Editable controls

| JSON key | Purpose |
|---|---|
| `mode` | `lid`, `thermal_lid`, or `natural` |
| `n` | Number of pressure/temperature cells per direction |
| `reynolds` | Active in lid and thermal_lid modes |
| `prandtl` | Thermal diffusivity relative to momentum diffusivity |
| `rayleigh` | Active in natural mode |
| `richardson` | Active in thermal_lid; set zero to disable buoyancy |
| `end_time` | Requested dimensionless simulation duration |
| `cfl` | Conservative explicit stability factor, maximum 0.5 |
| `dt_max` | Additional upper bound on the time step |
| `steady_tolerance` | Maximum u/v/theta change per time unit used for stopping |
| `max_steps` | Hard step cap; incomplete runs report failure |
| `samples` | Approximate number of saved history/animation samples |

Unused dimensionless groups remain in the common JSON schema but do not affect
the selected mode. Changing only `rayleigh` in lid mode, for example, does nothing.
For thermal_lid, Re, Pr, and Ri imply Ra = Ri Re² Pr, rather than using the separate
`rayleigh` entry. Keep parameters in modest laminar regimes for these examples.

## Numerical method

A Cartesian MAC layout puts pressure and temperature at cell centers, horizontal
velocity on vertical faces, and vertical velocity on horizontal faces. Normal
velocity is zero on cavity walls. Ghost values impose tangential no-slip and lid
motion; the moving-lid corner discontinuity remains a benchmark idealization.

Each step computes an explicit viscous/advective predictor, solves a Neumann pressure
Poisson problem with one reference value, and projects velocity onto the discrete
divergence-free space. The Poisson factorization is reused. Momentum convection uses
first-order upwind derivatives in advective form; it is not a kinetic-energy-
conserving discretization. Temperature uses conservative face fluxes. Time stepping
is forward Euler; diffusion is centered and second order in the interior. Numerical
upwind diffusion and projection boundary errors limit accuracy.

The adaptive step uses both convection and a conservative boundary-aware diffusion
bound. There is no nonlinear implicit solver. Doubling n increases both unknown
count and the number of explicit time steps; large grids can become expensive.

## Results

Each case writes:

| File | Contents |
|---|---|
| `plot.png` | Velocity, streamlines, pressure/temperature, profiles, convergence |
| `evolution.gif` | Startup animation when `--animate` is requested |
| `fields.csv` | Final cell-centered dimensionless fields |
| `fields.npz` | Final arrays including face velocities |
| `snapshots.npz` | Sampled u/v/temperature fields and times |
| `history.csv` | Divergence, kinetic energy, change rate, heat transfer |
| `result.json` | Inputs, termination status, metrics, environment, runtime |

`steady_tolerance_reached` means the field-change criterion held for 30 consecutive
steps. It does not mean mesh independence. `end_time_reached` means the requested
transient interval completed, without claiming steady state. `max_steps_reached`
saves an incomplete result and the command exits with an error.

Nusselt number measures wall heat transfer relative to pure conduction; Nu=1 is
the linear conduction limit. At steady state, hot-wall and cold-wall heat rates
should agree. During a transient their difference accounts for stored thermal
energy. `max_energy_step_defect` checks that discrete balance, not real-world
accuracy. Divergence is calculated from face velocities, not interpolated plot data.

## Experiments to try

1. Lid case: compare Re=50 and 100. Inspect circulation and centerline profiles.
2. Heated lid: compare Ri=0, 0.1, and 1 at fixed Re and Pr to isolate buoyancy.
3. Natural convection: compare Ra=1000 and 10000 at Pr=0.71.
4. Compare n=16, 32, and 48 before accepting quantitative results.
5. Halve cfl and dt_max to investigate time-step sensitivity at a fixed final time.

The supplied `--studies` command runs natural-convection refinement and a fixed-time
lid time-step comparison. See `VALIDATION.md` for measured results and limitations.

## References

- [FEATool developer tutorial: natural convection](https://www.featool.com/doc/Multiphysics_04_natural_convection1): benchmark setup and reported reference mean Nu=1.118 at Ra=1000, and 2.243 at Ra=10000, Pr=0.71.
- [De Vahl Davis, 1983](https://doi.org/10.1002/fld.1650030305): original square-cavity natural-convection benchmark study.
- [Ghia, Ghia & Shin, 1982](https://doi.org/10.1016/0021-9991%2882%2990058-4): a benchmark source for future detailed lid-cavity profile comparison; no claim of full Ghia validation is made in this pack.

The physical layout matches the natural-convection benchmark, but our solver uses
free-fall scaling and a lower-order advection scheme. Convert velocity/time scales
before comparing quantities reported using other nondimensionalizations.
