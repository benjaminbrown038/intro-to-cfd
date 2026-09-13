# Verification status

All five Python entry points were executed with the supplied default cases. Numerical verification is reproducible with `python3 tests/verify.py`; recorded results are in verification-results.json. The runtime versions are in environment.json.

## Checks passed

- Tier 1: analytical sine decay, boundary values, symmetry, positivity, explicit stability rejection, and approximately second-order combined refinement.
- Tier 2: periodic scalar mass conservation, decreasing analytical error for both advective schemes, negative velocity, pure advection, and pure diffusion.
- Flow solver: zero forcing/zero boundary velocity returns rest; channel velocity converges toward its analytical profile; cell continuity and linear residuals are small.
- 3D operator: manufactured divergence-free flow with three nonzero components and known body forcing. Relative cell-center velocity errors decreased from about 19.7% to 8.4% to 4.7% on 4, 6, and 8 cells per axis.
- Adaptation: inserted planes preserve existing edges, produce positive nonuniform widths, and maintain small cell divergence after solving.

The manufactured flow uses velocity = curl((1,2,3)*phi)*0.001, with phi=sin(pi*x)^2*sin(pi*y)^2*sin(pi*z)^2. Its no-slip boundaries and analytical force are evaluated in tests/verify.py. This checks the 3D operator independently of the adaptive cavity reference.

## Limits of the evidence

MATLAB/Octave is unavailable, so MATLAB code has not been executed or cross-checked numerically. Adaptive error plots use a finite numerical reference and interior sampling, not exact global error. The penalty-cylinder force has not been benchmarked. The examples are steady creeping flow; tests do not establish correctness for omitted inertial, turbulent, or transient physics.
