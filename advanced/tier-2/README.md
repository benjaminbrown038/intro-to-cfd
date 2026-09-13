# Tier 2: 1D advection–diffusion

## Implemented

Transport a periodic sine profile with constant speed and diffusivity. Compare upwind advection and Lax–Wendroff advection; both use explicit centered diffusion. The exact solution is a translated, exponentially decaying sine wave.

The conservative flux difference preserves total scalar mass. Positive and negative speeds, pure advection, and pure diffusion are supported. The selected time step satisfies a sufficient explicit stability bound. Lax–Wendroff is not positivity preserving in general; this smooth test is not a discontinuous-front benchmark. Diffusion is first-order in time; do not describe the combined scheme as uniformly second-order.

## Run

From the repository root, after installing requirements:

```bash
python3 advanced/tier-2/python/run_case.py
```

Optional flags: `--case PATH` and `--output PATH`. Default shared input: `cases/transport.json`.

For MATLAB, enter this tier's `matlab` directory and run `clear run_case` followed by `run_case`. The Python run was executed; the MATLAB implementation has not been executed here.

## Saved results

convergence.csv, per-scheme profile CSVs, solution.png, and run metadata.

Python saves under `results/python/`; MATLAB saves under `results/matlab/`. Fields are saved as NPZ in Python and MAT in MATLAB. Scalar comparisons use CSV. Metadata includes actual inputs; flow metrics distinguish assembly from sparse linear solve time.

## Experiments

Change one input at a time. Keep model assumptions valid and repeat numerical checks. Examine accuracy, conservation, and solve cost together. Higher mesh density alone does not establish correctness.
