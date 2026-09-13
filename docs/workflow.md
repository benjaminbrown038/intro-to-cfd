# Development workflow

1. Define one case, its assumptions, and a reference answer.
2. Implement a small runnable Python example.
3. Verify conservation, boundary conditions, and convergence.
4. Save plots and measurements under the tier's results folder.
5. Implement the equivalent MATLAB example using the same case.
6. Compare numerical results within a stated tolerance.
7. Update the tier README with actual run instructions and completion status.

## When the code grows

Within each language folder, introduce `meshing`, `physics_engine`, `solvers`, `postprocessing`, and `orchestration` modules as working code requires them. Avoid duplicating an entire platform structure before there is code to organize.

## Naming

Use lowercase snake_case for code files and functions; keep the tier folder names shown here. Use descriptive case names such as `laminar_channel`. Avoid ambiguous filenames such as `test_final2`.

## Dependencies

No numerical dependencies are required for this folder starter. Add Python requirements and MATLAB version/toolbox information when the first implementation exists.
