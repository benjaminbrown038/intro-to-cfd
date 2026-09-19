# Tier 3: 2D creeping-flow cavity

## Implemented

Solve steady momentum and continuity together on a staggered finite-volume grid. The upper wall moves with u=U*sin(pi*x)^2 on the unit square; other walls are stationary.

No advective acceleration or time integration is included. With the supplied density, viscosity, speed, and unit length, Re=0.01. The smooth lid avoids the corner jump of the classic constant-speed lid benchmark, so those benchmark values are not directly comparable.

## Model and boundary conditions

Steady incompressible Stokes flow without a body force:

$$
-\mu\nabla^2\mathbf{u}+\nabla p=\mathbf{0}, \qquad \nabla\cdot\mathbf{u}=0.
$$

For the square of side $L=1$ m, the upper wall moves smoothly:

$$
\mathbf{u}(x,L)=\left(U\sin^2\left(\frac{\pi x}{L}\right),0\right), \qquad \mathbf{u}=\mathbf{0}\ \text{on the other walls}.
$$

Density is used to characterize the creeping-flow regime:

$$
\mathrm{Re}=\frac{\rho UL}{\mu}=0.01\quad\text{for the supplied case}.
$$

Here $\mu$ is dynamic viscosity, $\rho$ is density, and $U$ is peak lid speed. Pressure needs a reference value; the solver pins one pressure cell.

## Run

From the repository root, after installing requirements:

```bash
python3 advanced/tier-3/python/run_case.py
```

Optional flags: `--case PATH` and `--output PATH`. Default shared input: `cases/cavity.json`.

For MATLAB, enter this tier's `matlab` directory and run `clear run_case` followed by `run_case`. The Python run was executed; the MATLAB implementation has not been executed here.

## Saved results

velocity/pressure fields, divergence, speed/quiver plot, and algebraic residual and timing metrics.

Python saves under `results/python/`; MATLAB saves under `results/matlab/`. Fields are saved as NPZ in Python and MAT in MATLAB. Scalar comparisons use CSV. Metadata includes actual inputs; flow metrics distinguish assembly from sparse linear solve time.

## Experiments

Change one input at a time. Keep model assumptions valid and repeat numerical checks. Examine accuracy, conservation, and solve cost together. Higher mesh density alone does not establish correctness.
