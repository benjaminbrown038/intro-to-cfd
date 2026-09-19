# Advanced CFD

A progression from scalar transport to coupled velocity–pressure solutions and
adaptive meshes in Python and MATLAB. These examples develop the numerical tools
needed to investigate fluid motion, geometry, conservation, and solution accuracy.

The first two tiers establish transport fundamentals. Tiers 3–5 solve **steady,
incompressible Stokes flow**, where viscous effects dominate and fluid inertia is
neglected. They do not yet solve general unsteady Navier–Stokes flow.

## Tier guide

| Tier | Implemented case | Main concepts | Questions to investigate |
|---|---|---|---|
| [1 — Diffusion](tier-1/README.md) | Transient 1D sine-profile decay | Explicit integration, spatial derivatives, analytical error | How do diffusivity, mesh spacing, and time step affect decay and accuracy? |
| [2 — Advection–diffusion](tier-2/README.md) | Periodic scalar transport | Conservative fluxes, upwind and Lax–Wendroff advection | How do transport speed and numerical diffusion affect the solution? |
| [3 — Cavity flow](tier-3/README.md) | 2D cavity with a smoothly moving lid | Staggered finite volumes, pressure–velocity coupling, continuity | How does wall motion produce circulation in creeping flow? |
| [4 — Channel and obstacle](tier-4/README.md) | Analytical channel check and circular resistance region | Geometry representation, penalty methods, profile verification | How does an approximate obstacle redirect the flow? |
| [5 — Adaptive meshes](tier-5/README.md) | 2D and 3D cavities on nonuniform tensor grids | Gradient indicators, directional refinement, matched-cost comparisons | Does concentrating cells improve the chosen velocity metric? |

## Folder organization

Each `tier-N` contains:

| Location | Purpose |
|---|---|
| `README.md` | Model assumptions, usage, outputs, and experiments |
| `cases/` | JSON inputs shared by the Python and MATLAB implementations |
| `python/run_case.py` | Python entry point |
| `matlab/run_case.m` | MATLAB entry point |
| `results/python/` | Generated Python results |
| `results/matlab/` | Generated MATLAB results |

Shared flow operators and drivers live in [common](../common/README.md), outside
this folder. Keep the repository structure intact so entry points can locate them.

## Run with Python

Open a terminal at the **repository root**, the directory containing `advanced/`,
`common/`, and `requirements.txt`. If you are currently inside `advanced/`, first
run `cd ..`.

The supplied dependency set targets Python 3.11 or newer.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

On Windows PowerShell, activate with `.venv\Scripts\Activate.ps1` instead.

Run all five supplied examples:

```bash
python run_all.py
```

Or choose a tier:

```bash
python advanced/tier-1/python/run_case.py
python advanced/tier-2/python/run_case.py
python advanced/tier-3/python/run_case.py
python advanced/tier-4/python/run_case.py
python advanced/tier-5/python/run_case.py
```

Each entry point supports `--case` and `--output`. For example:

```bash
python advanced/tier-4/python/run_case.py \
  --case advanced/tier-4/cases/channel.json \
  --output advanced/tier-4/results/channel_baseline
```

Copy a case JSON before editing it, and choose a distinct output directory for
comparisons. Rerunning with the same output location replaces named output files.

## Run with MATLAB

With MATLAB's current folder set to the repository root:

```matlab
run_all
```

For a single tier, change into its `matlab` directory:

```matlab
cd advanced/tier-3/matlab
clear run_case
run_case
```

All tiers intentionally use the same entry-point name. Avoid adding every tier's
MATLAB folder to the path at once. The provided implementations use base MATLAB
functions; their execution has not been verified in the supplied build environment.

## Understand the physics before changing inputs

### Tiers 1–2: transport building blocks

Tier 1 solves diffusion of a normalized scalar with fixed zero endpoint values.
Tier 2 adds constant-speed transport with periodic boundaries. Both have analytical
sine-wave solutions for comparison. A transported scalar is not automatically a
complete velocity, temperature, or concentration model: its interpretation requires
appropriate units, coefficients, and physical assumptions.

Tier 1 uses first-order time integration and second-order spatial differences.
Its default refinement changes time step with grid spacing, so the observed
convergence includes both temporal and spatial effects. Tier 2's combined
advection–diffusion method should not be described as uniformly second-order.

### Tiers 3–5: creeping flow

These tiers solve momentum and mass conservation together:

$$
-\mu\nabla^2\mathbf{u}+\nabla p=\mathbf{f}, \qquad \nabla\cdot\mathbf{u}=0.
$$

Viscosity is constant within a run. The absence of momentum advection and time
integration is a modeling assumption, not a switch enabled by raising the inlet
or lid speed. Changing density does not add inertial physics to a Stokes solver.
Check the Reynolds number and the case's assumptions before interpreting results.

Tier 3 uses a **smooth lid-speed distribution**. Its results are not directly
comparable to a standard cavity benchmark with a constant-speed lid and abrupt
corner velocity changes.

Tier 4 first checks a body-force-driven channel with matching parabolic velocities
prescribed at both ends. Its obstacle is a **circular resistance region**, not a
boundary-fitted impermeable cylinder. Residual flow inside the region is possible.
The reported penalty force is not a validated drag coefficient, and this case does
not predict vortex shedding.

Tier 5 inserts whole coordinate planes in selected intervals. This creates a
conforming nonuniform tensor mesh, not local octree or block refinement. Each new
mesh is solved as a fresh steady problem; no transient state is transferred.

## Equations by tier

### Tier 1



$$
\frac{\partial\phi}{\partial t}=\alpha\frac{\partial^2\phi}{\partial x^2}, \qquad 0\leq x\leq L.
$$

See [Tier 1](tier-1/README.md) for boundary conditions, numerical formulas, and limitations.

### Tier 2



$$
\frac{\partial\phi}{\partial t}+a\frac{\partial\phi}{\partial x}=\alpha\frac{\partial^2\phi}{\partial x^2}, \qquad \phi(x+L,t)=\phi(x,t).
$$

See [Tier 2](tier-2/README.md) for boundary conditions, numerical formulas, and limitations.

### Tier 3



$$
-\mu\nabla^2\mathbf{u}+\nabla p=\mathbf{0}, \qquad \nabla\cdot\mathbf{u}=0.
$$

See [Tier 3](tier-3/README.md) for boundary conditions, numerical formulas, and limitations.

### Tier 4



$$
-\mu\nabla^2\mathbf{u}+\nabla p=\mathbf{f}, \qquad \nabla\cdot\mathbf{u}=0.
$$

See [Tier 4](tier-4/README.md) for boundary conditions, numerical formulas, and limitations.

### Tier 5



$$
-\mu\nabla^2\mathbf{u}+\nabla p=\mathbf{0}, \qquad \nabla\cdot\mathbf{u}=0.
$$

See [Tier 5](tier-5/README.md) for boundary conditions, numerical formulas, and limitations.

## What to inspect in the results

Output files vary by tier. Common products include:

| Output | What to examine |
|---|---|
| Scalar or velocity profiles | Agreement with an analytical solution or expected behavior |
| Pressure and velocity fields | Circulation, flow diversion, and spatial gradients |
| CSV comparison tables | Error, conservation, mesh size, and cost |
| PNG plots | Profiles, mesh layout, speed fields, and refinement trends |
| NPZ or MAT field files | Numerical arrays for further analysis, where produced |
| JSON metadata and metrics | Actual case inputs, residuals, and recorded timings |

A small linear-system residual means the discrete equations were solved closely.
It does **not** establish that the mesh resolves the physical solution accurately.
Likewise, conservation and visually smooth plots do not replace a refinement study.

## Verification and mesh studies

Run the supplied Python numerical checks from the repository root:

```bash
python tests/verify.py
```

The supplied archive records execution of all five Python examples and successful
numerical checks; this README does not imply a fresh rerun of those simulations.
See [verification status](../docs/verification.md) for the recorded evidence and limits.

Use this workflow when extending a case:

1. Run the supplied baseline and save its inputs and outputs.
2. State the quantity of interest: profile error, flow rate, sampled velocity, or another defined metric.
3. Change one physical or numerical input at a time.
4. Compare at least three meshes while keeping the physical problem fixed.
5. For transient transport, check time-step sensitivity separately.
6. Examine conservation, reference error, and computational cost together.
7. Document remaining sensitivity and the limits of the model.

Tier 5 compares adaptive and uniform meshes with matching per-axis cell counts.
Its default reference is a finer **numerical** solution, not an exact answer.
The plotted velocity difference is evaluated on a fixed interior sample lattice;
it is not a whole-domain or wall-stress error estimate. Increase reference
resolution separately before making stronger accuracy claims.

## Relationship to the later case packages

The separately distributed `cfd-materials-cases` and `cfd-cavity-cases` packages
have their own entry points, dependencies, documentation, and tier labels. Their
numbering does not redefine the original `advanced/tier-1` through `tier-5` above.
The cavity mesh-study add-on belongs to `cfd-cavity-cases`.

Keep those packages intact and follow their own READMEs when installed alongside
this repository. Do not copy their commands into the original tier entry points.

## Future work

Potential extensions to this original tier sequence include:

- Time-dependent Navier–Stokes flow with momentum advection.
- Developing channel flow and boundary-fitted obstacle geometry.
- Coupled heat transfer and temperature-dependent fluid properties.
- Local adaptive refinement with conservative state transfer.
- Compressible gases, multiphase flow, and fluid–structure interaction.

These are development directions, not implemented capabilities of Tiers 1–5.
Each new model should include a focused benchmark and updated scope documentation.

## Further documentation

- [Repository overview](../README.md)
- [Numerical methods and assumptions](../docs/numerical-methods.md)
- [Verification status](../docs/verification.md)
- [Development workflow](../docs/workflow.md)
- [Complexity roadmap](../docs/complexity-roadmap.md)
