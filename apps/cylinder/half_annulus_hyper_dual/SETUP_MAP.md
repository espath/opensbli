# Half Annulus Setup Map

This file is specific to:
- `apps/cylinder/half_annulus_hyper_dual/half_annulus_hyper_dual.py`
- `apps/cylinder/half_annulus_hyper_dual/half_annulus_ns_baseline.py`

## Mesh File Expectations

Canonical mesh:
- `apps/cylinder/grid_generation/data_upper_half.h5`
- expected:
  - `shape = (251, 371)` with halo 5
  - `size = [241, 361]` (interior)

Interpretation:
- `dir0 = radial`
- `dir1 = angular`

## Boundary Mapping (Current Intent)

- `dir0, side0`: inner wall (isothermal no-slip)
- `dir0, side1`: outer farfield (Dirichlet freestream)
- `dir1, side0/1`: cut lines (`y=0`)

Physical preference for cut lines:
- symmetry BC on both sides

Current toolchain issue:
- direct `SymmetryBC` in this case generates invalid stencil names (compile failure).

## Current Working Setup

Validated working at `Ma=2.0`, `Re=2000`, `np=(241,361)`:

- Mesh:
  - regenerate with `beta_r=1.2` using `generate_cylinder_annulus.py`
  - this uses analytic halo construction in `(r,theta)` space.
- BCs:
  - `dir0,side0`: `IsothermalWallBC`
  - `dir0,side1`: `DirichletBC` (freestream)
  - `dir1,side0`: `SymmetryBC`
  - `dir1,side1`: `SymmetryBC`
- Time step:
  - `dt=5e-6` (default in `run_half_annulus_hyper_dual.sh`)

Observed:
- NS and hyper-dual (`Ml=0`) are both stable through 100 iterations.
- Interior non-finite count: `0`.
- NS/hyper parity at 100 iterations: ~`1e-15`.

## Debug Priorities

1. Keep this setup as the half-annulus reference baseline.
2. Port the same mesh/BC/dt discipline to quarter and full annulus.
3. Re-run smoke matrix after any mesh-generator change.
