# Working on the code

1. Run one tier with its supplied case before editing.
2. Change one physical or numerical input at a time.
3. Run `python3 tests/verify.py` after modifying numerical operators.
4. Inspect plots, conservation, and reference errors together.
5. Keep MATLAB and Python cases consistent; MATLAB parity still needs execution in MATLAB.
6. Update the tier README whenever model scope changes.

Shared flow equations live in common/python/stokes.py and common/matlab/cfd_stokes.m. Drivers, adaptation comparisons, and visualization live alongside those files. Each tier's entry point locates these dependencies automatically.

Python dependencies are pinned in requirements.txt; Python 3.11+ is required. Example output values and timings reflect the environment recorded in docs/environment.json. Results are generated under each tier and ignored by Git by default.
