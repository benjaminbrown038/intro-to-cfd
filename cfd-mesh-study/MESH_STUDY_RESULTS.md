# Completed mesh study

Nine runs completed: three cavity cases on 16×16, 32×32, and 48×48 grids.
All nine reached their steady tolerance. Thirteen tests passed: the original
nine solver tests and four new study/interpolation/assessment tests.

Changes from 32×32 to 48×48:

| Case | Velocity field change | Kinetic energy change | Hot-wall heat-transfer change |
|---|---:|---:|---:|
| tier_6_lid_cavity | 3.018% | 4.699% | — |
| tier_8_heated_lid | 2.460% | 3.475% | 0.132% |
| tier_9_natural_convection | 1.158% | 0.710% | 0.526% |

The 1% threshold is a study setting, not a universal acceptance criterion. These
results show why heat transfer alone is insufficient to judge every part of the
solution. The reports flag remaining changes; they do not claim exact error or
mesh independence. Further refinement should focus on the quantity you need.

The runs keep physics and cavity geometry fixed and use a common time-step cap
per case. Temporal sensitivity still needs a separate check before precise spatial
error claims. View each case's mesh_comparison.png for profiles and cost; use
mesh_summary.csv for exact values. Runtime is measured once and is sensitive to
machine load, so it is not a hardware benchmark.

## Source provenance

This additive extension was verified against the previously supplied cavity solver.
SHA-256 checksums of the base modules used for the sample runs:

- `solver.py`: `cba3047619693280f11a0c985728e7c434eb4793fa9a52538f1b8b54b98143ce`
- `__main__.py`: `80c681b20f771723d97951fcaa6982933b32e831a5348c2ba0314b07c6a40742`
