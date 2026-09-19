# Tier 2: 1D advection–diffusion

## Implemented

Transport a periodic sine profile with constant speed and diffusivity. Compare upwind advection and Lax–Wendroff advection; both use explicit centered diffusion. The exact solution is a translated, exponentially decaying sine wave.

The conservative flux difference preserves total scalar mass. Positive and negative speeds, pure advection, and pure diffusion are supported. The selected time step satisfies a sufficient explicit stability bound. Lax–Wendroff is not positivity preserving in general; this smooth test is not a discontinuous-front benchmark. Diffusion is first-order in time; do not describe the combined scheme as uniformly second-order.

## Model and discretization

Constant-speed transport with diffusion:

$$
\frac{\partial\phi}{\partial t}+a\frac{\partial\phi}{\partial x}=\alpha\frac{\partial^2\phi}{\partial x^2}, \qquad \phi(x+L,t)=\phi(x,t).
$$

The supplied sine profile has the exact solution:

$$
\phi(x,t)=1+\frac{1}{2}\exp\left[-\alpha\left(\frac{2\pi}{L}\right)^2t\right]\sin\left(\frac{2\pi(x-at)}{L}\right).
$$

The conservative update is:

$$
\phi_i^{n+1}=\phi_i^n-\frac{\Delta t}{\Delta x}\left(F_{i+1/2}^n-F_{i-1/2}^n\right).
$$

Upwind advection with centered diffusive flux uses:

$$
F_{i+1/2}^{\mathrm{UP}}=\max(a,0)\phi_i+\min(a,0)\phi_{i+1}-\alpha\frac{\phi_{i+1}-\phi_i}{\Delta x}.
$$

The comparison uses Lax–Wendroff advection with the same diffusion term:

$$
F_{i+1/2}^{\mathrm{LW}}=\frac{a}{2}(\phi_i+\phi_{i+1})-\frac{a^2\Delta t}{2\Delta x}(\phi_{i+1}-\phi_i)-\alpha\frac{\phi_{i+1}-\phi_i}{\Delta x}.
$$

All flux values above are evaluated at time level $n$. Define the signed Courant number and diffusion number:

$$
C=\frac{a\Delta t}{\Delta x}, \qquad D=\frac{\alpha\Delta t}{\Delta x^2}.
$$

The code selects a time step with the sufficient bound below (the configured safety factor further reduces it):

$$
|C|+2D\leq 1.
$$

Here $a$ is speed, $\alpha$ is diffusivity, and $L$ is the periodic domain length.

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
