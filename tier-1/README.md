# Tier 1: 1D diffusion

## Status

Transient sine-decay diffusion is implemented in Python and MATLAB. Python was executed and checked; MATLAB was not executed in the build environment. Steady diffusion and independent time-refinement exercises remain future additions. Tiers 2–5 now have runnable examples; see the repository README for their scope.

## Model

Equation: du/dt = alpha * d²u/dx², on 0 <= x <= L.
Boundary values: u(0,t) = u(L,t) = 0.
Initial condition: u(x,0) = sin(pi*x/L).
Analytical solution: u(x,t) = sin(pi*x/L)*exp(-alpha*(pi/L)^2*t).

The transported scalar is normalized and dimensionless. Length is in meters, time in seconds, and diffusivity in m²/s. This is a diffusion building block, not a full fluid velocity/pressure solver.

## Numerical method

Centered second-order spatial differences and forward Euler time integration. The explicit diffusion number alpha*dt/dx² must be at most 0.5. The code reduces dt as needed to reach the final time exactly.

The supplied refinement study uses dt proportional to dx², so both the spatial error and first-order temporal error scale with dx². Observed second-order convergence here is a combined refinement result; it does not imply second-order time integration.

## Run Python from the repository root

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
python3 advanced/tier-1/python/run_case.py
```

Optional: pass `--case PATH` for another case JSON, or `--output PATH` for a different output directory. Use increasing integer node counts with at least 3 nodes and positive physical parameters. The supplied case has 21, 41, 81, and 161 nodes.

## Run MATLAB from the repository root

```matlab
addpath('advanced/tier-1/matlab');
run_case
```

Uses base MATLAB functions including jsondecode, table, writetable, and plotting. No specialist toolbox is used. MATLAB execution is unverified here.

## Outputs

Both implementations save final numerical and analytical profiles and convergence tables under `results/python` or `results/matlab`. Each also saves solution and convergence PNG plots and run metadata. The Python outputs included in this ZIP are example results from the supplied case; rerunning replaces those named outputs.

## Verification performed

Python: zero boundary values, symmetry, nonnegative bounded values, stability parameter rejection, and analytical convergence. L2 errors were approximately 3.3132e-4, 8.2548e-5, 2.0620e-5, and 5.1538e-6. Observed refinement orders were 2.005, 2.001, and 2.000.

Integration timings measure only the time-stepping loop, excluding plotting, input/output, and setup. These vary by hardware and run.

## Experiments

1. Change diffusivity and observe how quickly the profile decays.
2. Change final time and compare numerical and analytical curves.
3. Compare diffusion numbers 0.1 and 0.4 at a fixed grid.
4. Try 0.6: the solver should reject the unstable choice.
5. Compare Python and MATLAB CSV values, excluding timing columns.
