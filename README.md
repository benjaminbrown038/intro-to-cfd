# intro-to-cfd

Progressive computational fluid dynamics in Python and MATLAB.

## What runs now

| Tier | Implemented example | Scope |
| --- | --- | --- |
| 1 | Transient 1D diffusion and analytical convergence | Explicit finite differences |
| 2 | Periodic advection–diffusion, upwind vs Lax–Wendroff advection | Conservative face fluxes; explicit diffusion |
| 3 | 2D cavity with a smoothly moving lid | Coupled steady Stokes velocity and pressure |
| 4 | Channel verification and a penalized circular obstacle | Steady creeping flow; approximate obstacle geometry |
| 5 | Adaptive 2D and 3D cavity simulations | Nonuniform tensor meshes; matched-cell-count uniform comparisons |

All five Python examples have been executed. Python numerical checks pass. MATLAB implementations are included for every tier, but have not been executed because MATLAB/Octave is unavailable in the build environment.

**Scope matters:** Tiers 3–5 solve steady incompressible Stokes flow, neglecting fluid inertia. They are not full time-dependent Navier–Stokes solvers. Tier 5 inserts whole coordinate planes based on computed velocity gradients; it does not implement local octree/block AMR. The roadmap retains these extensions for later work.

## Quick start: Python

Use Python 3.11 or newer. Open a terminal in this extracted repository folder:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
python3 run_all.py
```

Run one tier instead:

```bash
python3 advanced/tier-5/python/run_case.py
```

Run numerical verification:

```bash
python3 tests/verify.py
```

## Quick start: MATLAB

Open the repository root in MATLAB and run:

```matlab
run_all
```

For a single tier, navigate into its `matlab` folder, run `clear run_case`, then `run_case`. Do not add all five tier folders to the MATLAB path: they intentionally share the entry-point name. Base MATLAB functions are used; no specialist toolbox is required. MATLAB results are unverified here.

## Organization

- `advanced/tier-N/python/`: runnable entry point and tier-specific Python code.
- `advanced/tier-N/matlab/`: matching MATLAB entry point/code.
- `advanced/tier-N/cases/`: shared JSON physical and numerical inputs.
- `advanced/tier-N/results/`: generated plots, fields, and measurements.
- `common/python/` and `common/matlab/`: shared 2D/3D Stokes solver and drivers.
- `tests/verify.py`: Python numerical checks, including a manufactured 3D flow.
- `docs/`: equations, scope, verification results, and development roadmap.

Python outputs from the supplied cases are included in this archive. Rerunning replaces the named outputs. Generated results are ignored by Git by default; select small reference outputs deliberately if you want to track them. No GitHub repository has been changed by creating this archive.

The `advanced` parent preserves the structure established alongside intro-to-fea; the first two tiers are foundations.

Start with [Tier 1](advanced/tier-1/README.md), or inspect [Tier 5](advanced/tier-5/README.md). Read [the numerical methods](docs/numerical-methods.md) before changing physical parameters.
