# intro-to-cfd

Progressive computational fluid dynamics in Python and MATLAB—from fundamental transport equations to 3D adaptive mesh simulations.

## Current status

This starter contains folders and implementation guidance. Numerical solvers are not implemented yet. Folder creation does not indicate completion of a tier.

## Learning progression

| Tier | Topic | Status |
| --- | --- | --- |
| [tier-1](advanced/tier-1/README.md) | 1D diffusion | Planned |
| [tier-2](advanced/tier-2/README.md) | 1D advection–diffusion | Planned |
| [tier-3](advanced/tier-3/README.md) | 2D incompressible cavity flow | Planned |
| [tier-4](advanced/tier-4/README.md) | Channel flow and cylinder flow | Planned |
| [tier-5](advanced/tier-5/README.md) | Adaptive mesh CFD, progressing to 3D | Planned |

## Organization

- `advanced/tier-N/python/`: Python implementation for that tier.
- `advanced/tier-N/matlab/`: equivalent MATLAB implementation.
- `advanced/tier-N/cases/`: shared physical inputs, units, meshes, and boundary conditions.
- `advanced/tier-N/results/`: locally generated plots and measurements.
- `docs/`: workflow and verification guidance.

The `advanced` parent keeps the layout familiar alongside intro-to-fea; tiers 1–2 build the foundations.

## First complete flow example

The first complete flow milestone is Tier 4 laminar channel flow. Tiers 1–3 provide supporting exercises and numerical techniques. Implement the Python example first, verify it, and then match its inputs and outputs in MATLAB.

## Using this starter

Extract the ZIP. For a new repository, use the extracted `intro-to-cfd` folder as its root. If you already cloned an empty repository, copy this folder's contents into that clone. Review any existing files before replacing them.

Each folder contains a README so Git can track it. Add real source files as functionality is implemented; no empty solver files are supplied.

See [workflow](docs/workflow.md) for conventions and [verification](docs/verification.md) for completion checks.

## Increasing complexity

Each tier contains smaller steps so physics, geometry, and numerical difficulty can grow gradually. See the [complexity roadmap](docs/complexity-roadmap.md).
