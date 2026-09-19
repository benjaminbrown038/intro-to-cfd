# Mesh study extension

This adds a mesh-study command to the existing `cfd-cavity-cases` package. It requires
that package's `solver.py` and `__main__.py`; it is not a standalone solver.

## Install the additive ZIP

Extract `cfd-mesh-study-addon.zip` into your `intro-to-cfd` repository and merge the
matching directories. The archive contains only these additions beneath
`cfd-cavity-cases`:

- `cavity_cases/mesh_study.py`
- `tests/test_mesh_study.py`
- `MESH_STUDY.md`
- `MESH_STUDY_RESULTS.md`
- `mesh_study_verification.log`
- `mesh_results/` with completed sample runs, plots, tables, and reports

It contains no replacement solver, case JSON, requirements, or existing README.
If your file manager offers to replace an entire directory instead of merging,
copy the new files into their matching directories individually. Preserve any
local edits you have made. No new dependencies are needed.

## Run

From your existing `cfd-cavity-cases` folder with its environment active:

```bash
python -m cavity_cases.mesh_study --all
```

To study just natural convection:

```bash
python -m cavity_cases.mesh_study --case cases/tier-9/case.json
```

Choose your own meshes, output folder, comparison threshold, and time limit:

```bash
python -m cavity_cases.mesh_study --all --grids 16 32 64 --out mesh_64
python -m cavity_cases.mesh_study --case cases/tier-8/case.json --end-time 120 --threshold-percent 0.5
python -m unittest discover -s tests -v
```

At least three distinct grids with n>=8 are required. Grids are sorted before
running. Larger meshes may require longer runs or a higher max_steps in the input
JSON. The default --end-time is 80 dimensionless time units, overriding the case's
end_time for this steady study. It is an upper limit: the solver stops earlier when
its existing steady tolerance is satisfied. A max-step limit fails the command;
the incomplete run remains saved. Runs start from the original initial condition.

Results with the same case name in the same output folder are replaced on rerun.
Use --out to preserve an earlier experiment. All physics and geometry settings
come from the case JSON and stay fixed across its meshes. Only resolution,
case/output name, a common time-step cap, and the study time limit are changed.

## What gets compared?

| Grid | Pressure / temperature cells |
|---|---:|
| 16×16 | 256 |
| 32×32 | 1,024 |
| 48×48 | 2,304 |

These are uniform Cartesian cells. Face velocities occupy a separate staggered
layout: 2n(n+1) stored face values, including prescribed boundary values. The study
refines the same square geometry; it does not test mesh shapes, boundary-fitted
geometry, adaptive refinement, or wall-normal stretching.

Each successive mesh is compared using:

- **Horizontal and vertical centerline velocity profiles** to see flow-pattern changes.
- **Integrated kinetic energy** to compare overall flow strength.
- **Mean hot-wall Nusselt number** to compare thermal performance in heated cases.
- **Velocity field difference**: finer-grid velocities interpolated onto coarse cell
  centers, then a combined u/v relative L2 difference normalized by the mapped
  finer velocity norm.
- **Temperature field difference**: RMS theta difference as a percentage of the
  imposed hot-to-cold temperature range. This is not a percentage of absolute kelvin.
- **Final divergence and field-change rate** to check mass conservation and stopping.
- **Runtime, steps, and cell count** to show computation cost.

Coarse cell centers lie inside the finer grid's cell-center bounds; interpolation
never extrapolates. Linear interpolation is not a conservative restriction, and
these differences are not exact-solution errors. Pressure is saved in each run but
is not used as an accuracy acceptance metric because pressure and wall stresses
have not been independently validated in this solver.

Relative scalar change is `100*abs(new-old)/abs(new)`. A zero denominator gives an
undefined result unless both values are effectively zero. The first grid has no
preceding-grid change. Reported changes compare adjacent levels, not always the
coarsest and finest grid.

## Time-step and steady-state controls

The solver is explicit. Refining the mesh normally reduces its stable time step,
which can mix temporal and spatial effects. This study computes one conservative
cap from the finest grid and applies it to every grid. The adaptive solver may
still reduce it; mean_dt and dt_cap are recorded to make this visible.

The cap is `min(input dt_max, cfl/(8 D n_fine² + 2 n_fine))`, where D is the largest
active dimensionless diffusivity. This reuses the solver's conservative diffusion
bound and provides a reference convective allowance. The adaptive stability bound
remains active at every step. A common cap reduces temporal contamination but is
not a time-step independence test.

For a separate temporal check, copy the JSON and set its dt_max below the reported
cap, then repeat the same meshes in a separate output directory. Halving the time
step should change the quantity of interest much less than the spatial refinement
changes. For a true transient study, compare solutions at the same physical time;
this extension is designed for steady-state comparisons and permits early steady
termination at different times.

## Reading the assessment

The default threshold is 1%. It checks the **finest pair** for kinetic energy and
velocity-field changes, plus Nu and temperature change in thermal cases. All grid
runs must also reach the steady criterion. The JSON records the steady-state and
change-threshold flags separately.

- **Above threshold:** refine further for the selected metrics.
- **Within threshold and all steady:** the finest pair agrees to the chosen level.
- **Any nonsteady run:** the steady comparison is inconclusive; extend the time limit.

Even a passing threshold does not prove mesh independence. Check another level,
benchmark quantities, temporal sensitivity, and whether the important physical
features are resolved. The solver has first-order upwind advection, so do not assume
second-order overall accuracy. This extension intentionally does not report a Grid
Convergence Index or a Richardson-extrapolated exact answer: those require additional
assumptions about systematic refinement and asymptotic behavior.

## Output files

Each case gets its own folder in `mesh_results`:

| File | Contents |
|---|---|
| `mesh_comparison.png` | Six panels: profiles, heat transfer/energy, changes, cost, diagnostics |
| `mesh_summary.csv` | One row per mesh with inputs, metrics, and differences |
| `assessment.json` | Threshold and steady-state flags with interpretation |
| `REPORT.md` | Readable per-case results |
| `runs/grid_N/` | Full per-grid fields, history, plot, and metadata |

`mesh_results/all_cases_mesh_summary.csv` combines all selected cases. All dimensional
interpretation uses the scaling documented in the original package README.

To keep regenerable outputs out of git, add `mesh_results/` (and other custom output
folders you use) to your local `.gitignore`. Keep the study code and documentation
under version control. The original solver's assumptions and limitations apply.
