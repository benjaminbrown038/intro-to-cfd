# Tier 5: Adaptive 2D and 3D creeping flow

## Implemented

Solve a smooth-lid cavity, compute velocity-gradient activity, mark coordinate intervals, insert their midpoint planes, and solve again. Repeat in both 2D and 3D. Each cycle also solves a uniform mesh with the same per-axis cell counts.

Refinement is anisotropic tensor-grid refinement: a selected interval is split across the entire domain. There are no hanging faces or octrees. Each mesh is solved from scratch as a steady problem, so there is no transient state transfer. Conservation is enforced using shared face fluxes on each mesh. The gradient indicator is heuristic, not a certified error bound.

## Model and mesh equations

Each 2D or 3D mesh solves steady incompressible Stokes flow:

$$
-\mu\nabla^2\mathbf{u}+\nabla p=\mathbf{0}, \qquad \nabla\cdot\mathbf{u}=0.
$$

For a control volume $P$, discrete continuity balances outward face fluxes:

$$
\sum_{f\in\partial P}(\mathbf{u}_f\cdot\mathbf{n}_f)A_f=0.
$$

Here $A_f$ is face area in 3D (edge length per unit depth in 2D), and $\mathbf{n}_f$ is the outward unit normal. The velocity-gradient indicator ranks coordinate intervals. Each marked interval gains a midpoint:

$$
x_{i+1/2}^{\mathrm{new}}=\frac{x_i+x_{i+1}}{2}.
$$

The same operation applies along the other coordinate directions; each inserted coordinate creates an entire mesh plane. Tensor-grid cell counts are:

$$
N_{\mathrm{cells}}=N_xN_y\quad\text{in 2D}, \qquad N_{\mathrm{cells}}=N_xN_yN_z\quad\text{in 3D}.
$$

The reported sampled velocity difference is:

$$
E_h=\frac{\left[\sum_j\left\|\mathbf{u}_h(\mathbf{x}_j)-\mathbf{u}_{\mathrm{ref}}(\mathbf{x}_j)\right\|_2^2\right]^{1/2}}{\left[\sum_j\left\|\mathbf{u}_{\mathrm{ref}}(\mathbf{x}_j)\right\|_2^2\right]^{1/2}}.
$$

Both fields are linearly interpolated onto the same fixed interior sampling lattice. The reference is numerical, so this metric includes reference and interpolation error.

## Run

From the repository root, after installing requirements:

```bash
python3 advanced/tier-5/python/run_case.py
```

Optional flags: `--case PATH` and `--output PATH`. Default shared input: `cases/adaptive.json`.

For MATLAB, enter this tier's `matlab` directory and run `clear run_case` followed by `run_case`. The Python run was executed; the MATLAB implementation has not been executed here.

## Saved results

all adaptive/uniform/reference fields, nonuniform 3D mesh plots, speed slices, comparison.csv, and adaptive_comparison.png.

Python saves under `results/python/`; MATLAB saves under `results/matlab/`. Fields are saved as NPZ in Python and MAT in MATLAB. Scalar comparisons use CSV. Metadata includes actual inputs; flow metrics distinguish assembly from sparse linear solve time.

## Experiments

Change one input at a time. Keep model assumptions valid and repeat numerical checks. Examine accuracy, conservation, and solve cost together. Higher mesh density alone does not establish correctness.

## Meaning of the error plot

The reference is a **13-cell-per-axis numerical solution**, not an exact solution. Error is the relative Euclidean norm of velocity differences at a fixed interior sampling lattice (9 points per axis from 0.15 to 0.85 in the unit domain), after linear interpolation. It does not measure whole-domain error or wall stress error. Reference, interpolation, and coarse-grid errors all contribute.

The default 3D adaptive sequence uses 125, 343, and 729 cells. In the executed run, sample differences were approximately 41.9%, 14.7%, and 7.5%. Uniform meshes with matching final cell counts gave approximately 41.9%, 22.5%, and 10.6%. These results demonstrate this case only; adaptation need not outperform uniform refinement for every case or indicator.

Increase `reference_cells` independently to assess reference sensitivity before drawing stronger conclusions. The sparse direct solver's memory demand rises quickly in 3D; the defaults are deliberately small.
