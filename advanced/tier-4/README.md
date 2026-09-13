# Tier 4: Channel verification and penalized cylinder

## Implemented

Compare channel velocity against u(y)=4*Umax*y*(1-y) for unit height. A streamwise body force 8*mu*Umax drives the flow, and matching parabolic velocities are prescribed at both ends. Then add a circular penalty region at the domain center.

The circular region is represented by a Brinkman resistance on a Cartesian mesh, not a fitted impermeable boundary. Penalty force is an approximate volume-integrated resistance, not a validated drag coefficient. Residual flow inside the region is expected. No developing-flow or unsteady vortex-shedding model is implemented.

## Run

From the repository root, after installing requirements:

```bash
python3 advanced/tier-4/python/run_case.py
```

Optional flags: `--case PATH` and `--output PATH`. Default shared input: `cases/channel.json`.

For MATLAB, enter this tier's `matlab` directory and run `clear run_case` followed by `run_case`. The Python run was executed; the MATLAB implementation has not been executed here.

## Saved results

channel profiles, comparison.csv, velocity/pressure fields, obstacle flow plot, and penalty-force estimates.

Python saves under `results/python/`; MATLAB saves under `results/matlab/`. Fields are saved as NPZ in Python and MAT in MATLAB. Scalar comparisons use CSV. Metadata includes actual inputs; flow metrics distinguish assembly from sparse linear solve time.

## Experiments

Change one input at a time. Keep model assumptions valid and repeat numerical checks. Examine accuracy, conservation, and solve cost together. Higher mesh density alone does not establish correctness.
