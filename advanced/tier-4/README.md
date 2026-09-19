# Tier 4: Channel verification and penalized cylinder

## Implemented

Compare channel velocity against u(y)=4*Umax*y*(1-y) for unit height. A streamwise body force 8*mu*Umax drives the flow, and matching parabolic velocities are prescribed at both ends. Then add a circular penalty region at the domain center.

The circular region is represented by a Brinkman resistance on a Cartesian mesh, not a fitted impermeable boundary. Penalty force is an approximate volume-integrated resistance, not a validated drag coefficient. Residual flow inside the region is expected. No developing-flow or unsteady vortex-shedding model is implemented.

## Model and geometry

The baseline channel solves:

$$
-\mu\nabla^2\mathbf{u}+\nabla p=\mathbf{f}, \qquad \nabla\cdot\mathbf{u}=0.
$$

For channel height $H=1$ m, the analytical profile and applied body force are:

$$
u(y)=4U_{\max}\frac{y}{H}\left(1-\frac{y}{H}\right), \qquad v=0, \qquad f_x=\frac{8\mu U_{\max}}{H^2}, \qquad f_y=0.
$$

The baseline has zero streamwise pressure gradient, no-slip channel walls, and matching parabolic velocities at both ends. Its flow rate per unit out-of-plane width is:

$$
Q^{\prime}=\int_0^H u(y)\,dy=\frac{2}{3}U_{\max}H.
$$

The obstacle case adds resistance inside the circular region:

$$
-\mu\nabla^2\mathbf{u}+\nabla p+\lambda\chi\mathbf{u}=\mathbf{f}, \qquad \nabla\cdot\mathbf{u}=0.
$$



$$
\chi(x,y)=\begin{cases}1,&(x-x_c)^2+(y-y_c)^2\leq R^2,\\0,&\text{otherwise}.\end{cases}
$$

Here $R$ and $(x_c,y_c)$ describe the circle. The coefficient $\lambda$ has units Pa·s/m²; finite resistance permits residual interior flow. A continuous analogue of the reported force estimate per unit depth is:

$$
\mathbf{F}^{\prime}_{\mathrm{penalty}}=\int_{\Omega}\lambda\chi\mathbf{u}\,dA.
$$

This is the resistance-force estimate on the obstacle; the corresponding force on the fluid has the opposite sign. It is not validated drag.

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
