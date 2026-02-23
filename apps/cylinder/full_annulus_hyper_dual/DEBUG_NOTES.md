# Full Annulus Debug Notes

## Hard Rule (from user)

- Do **not** spend time checking dual-vs-NS parity for annulus debugging.
- Treat hyper-dual equations as validated (airfoil parity already established).
- Focus debugging on:
  1) mesh generation/layout/metadata consistency,
  2) boundary-condition direction/side mapping,
  3) topology/seam treatment.

## Persistent References

- General setup reference:
  - `apps/cylinder/GEOMETRY_TOPOLOGY_BC_PLAYBOOK.md`
- Half-annulus case map:
  - `apps/cylinder/half_annulus_hyper_dual/SETUP_MAP.md`

## Objective
Stabilize `full_annulus_hyper_dual` for 100 iterations without NaNs in physical domain.

## Change Log

1. Added WENO metric-class ordering fix in solver setup.
2. Added stale output cleanup in run script (`rm -f opensbli_output*.h5`).
3. Added full WENO periodic exchanges (`[-5,5]`) and explicit full-swap periodic kernels.
4. Added optional SFD + binomial filters to match cylinder stabilization pattern.
5. Mesh generator updates:
   - Full-circle angle endpoint handling (`2*pi` endpoint excluded).
   - Periodic halo wrapping for full-circle theta direction.
   - Optional transpose output.
   - Transpose-aware `npoints` metadata fix.
6. Run script now auto-detects default `NP0/NP1` from `data.h5` shape (`shape[1], shape[0]`) to avoid manual mismatch.
7. Solver script now reads mesh shape from `data.h5` and defaults `block0np0/block0np1` accordingly (same style as `supersonic_cylinder.py`).
8. Generator was re-tested with transposed output (`shape=(731,251)`) and corrected transpose metadata (`size=[721,241]`).
9. Run script now auto-derives `NP0/NP1` from `data.h5` shape when not explicitly provided.
10. Added missing `conservative=conservative` in `SimulationBlock(...)` to match working cylinder scripts.
11. Added `shock_factor` simulation constant for parity with cylinder setup.
12. Updated `NP0/NP1` auto-detection to prefer HDF5 `size` metadata over raw dataset shape.
13. Audited non-cylinder curvilinear examples (`airfoil_MB_2D_hyper_dual.py`, `transonic_MB.py`, `euler_wave_curvilinear.py`, `turbulent_counter_flow.py`) and added missing metric-velocity substitution:
    - `Eq(U_i, D_i_j*u_j)` in optional substitutions.

## Test Outcomes (recent)

- `ENABLE_WENO=1 ENABLE_SFD_FILTER=1 ENABLE_BINOMIAL_FILTER=1 NITER=100`:
  - Still fails NaN-check at iteration 100 in halo indices (e.g. `(-4,-5)`, `(-2,-5)`).
- With shape-derived `NP0/NP1=(251,731)` (from current data shape), still fails at iter 100 (`(-4,-5)`).
- With metadata-derived `NP0/NP1=(721,241)`, still fails at iter 100 (`(-4,-5)`).
- With NaN-check disabled and save at iteration 100:
  - `opensbli_output_000100.h5` still contains large NaN regions in interior.
- NS smoke variant (`full_annulus_ns_smoke.py`) showed similar behavior, indicating issue is not hyper-dual-only.

## Additional Findings (current session)

1. Cross-mesh smoke using the same `full_annulus_hyper_dual` driver:
   - `apps/cylinder/incompressible_cylinder/data.h5`: stable for 100 iterations when `PERIODIC_DIR=1`.
   - `apps/cylinder/grid_generation/data_full.h5`: NaN by iteration 100.
2. Filter isolation on annulus mesh:
   - Fails even with `ENABLE_WENO=0 ENABLE_SFD_FILTER=0 ENABLE_BINOMIAL_FILTER=0`.
3. NaN-check behavior:
   - With `ENABLE_NAN_CHECK=1`, first reported NaNs appear in halo indices (e.g. `(-2,-2)`).
   - With `ENABLE_NAN_CHECK=0`, run completes 100 iterations, but output still develops widespread NaNs (not only halos), so this is not just a diagnostics false positive.
4. Mesh-layout check:
   - `data_full.h5` (transposed) uses `shape=(731,251), size=[721,241]`.
   - `data_full_notrans.h5` uses `shape=(251,731), size=[241,721]`.
   - Both layouts still destabilize by iteration 100 under tested mappings.
5. Tried stabilizing periodic depth handling:
   - `PeriodicBC(..., full_depth=True)` on main periodic boundaries and WENO periodic swaps.
   - Improves initial finite ratio slightly in some runs but does not prevent eventual collapse.
6. Tried alternative wall model:
   - `IsothermalWall_ZeroPressureGradBC` at inner wall.
   - Did not resolve instability for the annulus meshes.
7. Half-annulus focused mesh/BC debug:
   - The active `data_upper_half.h5` had been in a transposed layout at one point (`shape=(371,251)`), which caused BC-direction confusion.
   - Regenerated clean non-transposed half mesh (`shape=(251,371), size=[241,361]`) and reran.
   - Corrected BC-direction mapping in half-annulus scripts:
     - dir0: radial (wall/farfield)
     - dir1: angular (cut lines)
   - This reduced early non-finite footprint significantly (from ~1.6e4 to ~1.8e3 cells at niter=1), but did not eliminate it.
8. Symmetry BC attempt on half-domain cut lines:
   - Root cause found in stencil name generation (`opensbli/core/kernel.py`), where negative stencil bounds produced invalid C identifiers (e.g. `stencil_..._2-1_...`).
   - Patched name encoding to C-safe signed tokens (`m#`, `p#`, `0`), so symmetry BC now compiles.
   - Post-fix half-annulus run (`niter=1`) with symmetry BC compiles and runs, but still shows non-finite interior values; setup is improved but not yet stable.
9. Mesh halo fill sensitivity:
   - Testing a copy-halo mesh variant worsened stability dramatically.
   - Original linear-extrapolated halo coordinates are currently less bad than copy-halo for this case.
10. Mesh-generator fix:
   - Updated annulus generator to build halo coordinates analytically in `(r, theta)` space (instead of linear extrapolation in `x/y`).
   - This reduced early non-finite footprint in half-annulus runs.
11. Stability envelope found for half-annulus:
   - Using `beta_r=1.2` mesh grading and `dt=5e-6`:
     - NS run stable to 100 iterations with zero non-finite interior cells.
     - Hyper-dual (`Ml=0`) run stable to 100 iterations with zero non-finite interior cells.
     - NS vs hyper-dual parity at 100 iterations: max diff ~ `7.8e-16` on interior fields.

## Working Configuration (Half Annulus)

- Mesh generator:
  - `generate_cylinder_annulus.py` with `--beta-r 1.2`
  - upper-half domain (`theta0=0`, `theta1=pi`)
- BCs:
  - `dir0`: inner wall (`IsothermalWallBC`), outer farfield (`DirichletBC`)
  - `dir1`: symmetry on both cut lines (`SymmetryBC`)
- Time step:
  - `dt=5e-6` for the tested `Ma=2.0, Re=2000, np=(241,361)`

## Current Working Hypothesis

The annulus case instability is tied to a setup mismatch between this full-ring topology and current BC/metric treatment (not to dual-velocity equations alone). The same solver stack is stable on at least one established curvilinear mesh (`incompressible_cylinder`) with matched direction mapping.

## Next Candidates

1. Re-check metric transformation inputs against direction mapping on this mesh (`Delta0/Delta1`, BC direction, periodic direction) using one consistent orientation only.
2. Build minimal no-filter/no-SFD/no-binomial run on same mesh to isolate first instability source.
3. If needed, temporarily reduce to half/quarter annulus where inflow/outflow placement is less ambiguous.
