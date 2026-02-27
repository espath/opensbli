# Quarter Annulus Debug Plan (Strict, Step-by-Step)

This file is the authoritative plan for debugging quarter-annulus instability.
No ad-hoc BC/geometry experiments outside this plan.

## Freeze Baseline (Step 1)

- Date: 2026-02-24
- Branch: `main`
- Commit: `dc170a4`
- Case file: `apps/cylinder/quarter_annulus_hyper_dual/quarter_annulus_hyper_dual.py`
- Runner: `apps/cylinder/quarter_annulus_hyper_dual/run_quarter_annulus_hyper_dual.sh`
- Mesh generator: `apps/cylinder/grid_generation/generate_cylinder_annulus.py`
- Mesh used: quarter annulus, `theta in [pi/2, pi]`, `r in [0.5, 20.0]`
- Resolution: `np=(241,241)`, halo `5`, `beta-r=1.2`
- Runtime params: `Ma=2.0`, `Ml=0.0`, `Re=1.0e4`, `dt=1.0e-5`, `niter=4000`, `save_every=1000`
- BC map:
  - `dir0,side0` (`theta=pi/2`): `PressureOutletBC`
  - `dir0,side1` (`theta=pi`): `SymmetryBC` (symmetry/slip)
  - `dir1,side0` (inner radius): `IsothermalWallBC`
  - `dir1,side1` (outer radius): `DirichletBC` (freestream)
- Known failure signature:
  - NaN at `index:(-2,-2)` around iteration ~3600
  - Same corner signature across multiple geometry/BC variants

## Objective

Fix reproducible instability in quarter-annulus run without changing physics model.

## Controlled Procedure

### Step 2: Reproduce Control

1. Clean quarter case outputs and generated wrappers.
2. Regenerate quarter mesh with fixed parameters.
3. Run control (`niter=4000`, `dt=1.0e-5`).
4. Confirm failure signature is unchanged.

If not reproducible, stop and re-freeze baseline.

### Step 3: Measure BC Write Overlap (No physics changes)

1. Inspect generated `opensbli_ops.cpp` BC kernel order.
2. Inspect BC kernel ranges for:
   - outlet (`dir0,side0`)
   - symmetry (`dir0,side1`)
   - wall (`dir1,side0`)
3. Build explicit overlap table for corner halos.

Expected: overlap at wall/symmetry corner and likely outlet touching corner halos.

### Step 4: Isolated Change A (Outlet Corner Trim)

Patch only `PressureOutletBC` application ranges to avoid writing corner halos.

- No changes to symmetry, wall, mesh, dt.
- Re-run control and compare:
  - fail iteration
  - fail index
  - early boundary behavior

### Step 5: If A fails, isolated change B (Corner ownership)

Revert A.
Apply explicit corner ownership rule (wall/symmetry corner owned by wall+symmetry, not outlet).

Re-run control and compare same metrics.

### Step 6: If B fails, isolated change C (corner compatibility update)

Revert B.
Add explicit corner compatibility update only at shared corner points.

Re-run control and compare.

## Acceptance Criteria

1. Survive beyond 4000 iterations under control settings, or
2. Same run no longer fails at `(-2,-2)` and shows materially delayed instability.

Must be reproducible in two consecutive runs.

## Guardrails

- Do not modify `hyper_dual_velocity_physics.py` during this debug.
- Do not change geometry topology while working this plan.
- Do not run multiple untracked hacks simultaneously.
- One code change at a time, always against control.

## Change Tracking (for easy restore)

- File: `opensbli/core/boundary_conditions/pressure_outlet.py`
  - Change: optional `include_corner_halos` toggle for `PressureOutletBC`.
  - Restore command: `git checkout -- opensbli/core/boundary_conditions/pressure_outlet.py`
- File: `opensbli/core/boundary_conditions/inviscid_wall.py`
  - Change: optional `include_corner_halos` toggle for `InviscidWallBC`.
  - Restore command: `git checkout -- opensbli/core/boundary_conditions/inviscid_wall.py`
- File: `apps/cylinder/quarter_annulus_hyper_dual/quarter_annulus_hyper_dual.py`
  - Change: BC corner toggles via env (`OUTLET_CORNER_HALOS`, `SYMMETRY_CORNER_HALOS`, `WALL_CORNER_HALOS`).
  - Change: optional corner debug monitor (`CORNER_DEBUG_MONITOR=1`) writing `quarter_annulus_corner_debug.log`.
  - Restore command: `git checkout -- apps/cylinder/quarter_annulus_hyper_dual/quarter_annulus_hyper_dual.py`

## Log Template (fill per run)

- Change ID:
- Files touched:
- Control or variant:
- Run command:
- Outcome:
  - Max iteration reached:
  - NaN index (if any):
  - First visible oscillation location:
  - Notes:

## Debug Log

- Change ID: `S2-control-repro`
- Files touched: none
- Control or variant: control
- Run command: `PYTHON_BIN=/Users/pmzle/Documents/dev/opensbli/espath312/bin/python3 NITER_VALUE=4000 SAVE_EVERY_VALUE=1000 DT_VALUE=1.0e-5 ./run_quarter_annulus_hyper_dual.sh`
- Outcome:
  - Max iteration reached: `3600` (NaN check at this print interval)
  - NaN index (if any): `(-2,-2)`
  - First visible oscillation location: not measured in this run
  - Notes: baseline failure signature reproduced after full clean + mesh regen.

- Change ID: `S3-bc-overlap-map`
- Files touched: none (inspection only)
- Control or variant: static inspection of generated code
- Run command: `rg/sed on opensbli_ops.cpp and opensbliblock00_kernels.h`
- Outcome:
  - Max iteration reached: n/a
  - NaN index (if any): n/a
  - First visible oscillation location: n/a
  - Notes:
    - BC call order per RK stage: `PressureOutlet(dir0,side0)` -> `Symmetry(dir0,side1)` -> `IsothermalWall(dir1,side0)` -> `Dirichlet(dir1,side1)`.
    - `PressureOutlet` iteration range: `i=0`, `j=[-2, np1+1]`; kernel writes `i={0,-1,-2}`.
    - `IsothermalWall` iteration range: `j=0`, `i=[-2, np0+1]`; kernel writes `j={0,-1,-2}`.
    - Therefore both BCs write lower-left corner halos, including `(-2,-2)`, in each stage.

- Change ID: `S4-outlet-corner-trim`
- Files touched: `opensbli/core/boundary_conditions/pressure_outlet.py`, `apps/cylinder/quarter_annulus_hyper_dual/quarter_annulus_hyper_dual.py`
- Control or variant: isolated variant (outlet only)
- Run command: `PYTHON_BIN=/Users/pmzle/Documents/dev/opensbli/espath312/bin/python3 OUTLET_CORNER_HALOS=0 NITER_VALUE=4000 SAVE_EVERY_VALUE=1000 DT_VALUE=1.0e-5 ./run_quarter_annulus_hyper_dual.sh`
- Outcome:
  - Max iteration reached: `3600`
  - NaN index (if any): `(-2,-2)`
  - First visible oscillation location: not measured in this run
  - Notes:
    - Generated outlet range changed as intended: `iteration_range_26_block0 = {0, 1, 0, block0np1}` (no corner halos in orthogonal direction).
    - Failure signature unchanged from control.

- Change ID: `S5-wall-corner-trim`
- Files touched: `apps/cylinder/quarter_annulus_hyper_dual/quarter_annulus_hyper_dual.py`
- Control or variant: isolated variant (wall corners only)
- Run command: `PYTHON_BIN=/Users/pmzle/Documents/dev/opensbli/espath312/bin/python3 OUTLET_CORNER_HALOS=1 WALL_CORNER_HALOS=0 NITER_VALUE=4000 SAVE_EVERY_VALUE=1000 DT_VALUE=1.0e-5 ./run_quarter_annulus_hyper_dual.sh`
- Outcome:
  - Max iteration reached: `3600`
  - NaN index (if any): `(-2,-2)`
  - First visible oscillation location: not measured in this run
  - Notes:
    - Generated wall range changed as intended: `iteration_range_28_block0 = {0, block0np0, 0, 1}` (no corner halos in orthogonal direction).
    - Failure signature unchanged from control.

- Change ID: `S6-dt-half-cfl-check`
- Files touched: none (runtime parameter sweep)
- Control or variant: timestep sensitivity
- Run command: `PYTHON_BIN=/Users/pmzle/Documents/dev/opensbli/espath312/bin/python3 NITER_VALUE=8000 SAVE_EVERY_VALUE=2000 DT_VALUE=5.0e-6 OUTLET_CORNER_HALOS=1 WALL_CORNER_HALOS=1 ./run_quarter_annulus_hyper_dual.sh`
- Outcome:
  - Max iteration reached: `7100`
  - NaN index (if any): `(-2,-2)`
  - First visible oscillation location: not measured in this run
  - Notes:
    - At `dt=5e-6`, failure occurs at simulation time `0.03550`.
    - Baseline (`dt=1e-5`) failed at simulation time `0.03600`.
    - Conclusion: instability is tied to physical-time evolution, not simply CFL from too-large `dt`.

- Change ID: `S7-corner-diagnostics`
- Files touched: `apps/cylinder/quarter_annulus_hyper_dual/quarter_annulus_hyper_dual.py`
- Control or variant: baseline with per-corner diagnostic monitor
- Run command: `... CORNER_DEBUG_MONITOR=1 MONITOR_EVERY_VALUE=1 ...`
- Outcome:
  - Max iteration reached: `3600` (NaN check trigger)
  - NaN index (if any): `(-2,-2)`
  - First visible oscillation location: corner neighborhood (tracked points)
  - Notes:
    - Diagnostic file: `block0_quarter_annulus_corner_debug.log`.
    - First non-finite monitored value appears at iteration `3529`, simulation time `0.03529`.
    - All monitored corner-neighborhood fields (`rho, rhou0, rhou1, rhoE, p, T` at points `(-2,-2),(-1,-1),(0,0),(0,1),(1,0),(1,1)`) become `NaN` in the same iteration.
    - At iteration `3528`, all monitored values are finite and smooth.

- Change ID: `S8-symmetry-corner-off`
- Files touched: none (runtime toggles only)
- Control or variant: disable symmetry corner halos
- Run command: `... OUTLET_CORNER_HALOS=1 SYMMETRY_CORNER_HALOS=0 WALL_CORNER_HALOS=1 ...`
- Outcome:
  - Max iteration reached: `3600`
  - NaN index (if any): `(-2,-2)`
  - First visible oscillation location: not measured in this run
  - Notes: no change from baseline failure signature.

- Change ID: `S9-symmetry-corner-owner`
- Files touched: none (runtime toggles only)
- Control or variant: corner ownership by symmetry only
- Run command: `... OUTLET_CORNER_HALOS=0 SYMMETRY_CORNER_HALOS=1 WALL_CORNER_HALOS=0 ...`
- Outcome:
  - Max iteration reached: `3600`
  - NaN index (if any): `(0,-2)`
  - First visible oscillation location: wall-halo line near outlet-side corner
  - Notes:
    - Failure location moved from `(-2,-2)` to `(0,-2)` under this ownership policy.
    - Indicates boundary-corner/near-corner update path affects where blow-up starts, but does not remove instability.

- Change ID: `S10-adiabatic-carpenter-wall`
- Files touched: `apps/cylinder/quarter_annulus_hyper_dual/quarter_annulus_hyper_dual.py`
- Control or variant: inner wall switched to `AdiabaticWall_CarpenterBC`
- Run command: `... WALL_BC_TYPE=adiabatic_carpenter NITER_VALUE=4000 DT_VALUE=1.0e-5 ...`
- Outcome:
  - Max iteration reached: `3300`
  - NaN index (if any): `(-2,0)`
  - First visible oscillation location: shifted to wall-adjacent halo line
  - Notes:
    - This is earlier than baseline fail (`~3600`).
    - Instability is not removed; it shifts and worsens in iteration count.

- Change ID: `S11-theta-pi-extrapolation-test`
- Files touched: `apps/cylinder/quarter_annulus_hyper_dual/quarter_annulus_hyper_dual.py`
- Control or variant: replaced `theta=pi` boundary with `ExtrapolationBC(order=0)`
- Run command: `THETA_PI_BC_TYPE=extrapolation NITER_VALUE=4000 SAVE_EVERY_VALUE=1000 DT_VALUE=1.0e-5 OUTLET_CORNER_HALOS=1 SYMMETRY_CORNER_HALOS=1 WALL_CORNER_HALOS=1 WALL_BC_TYPE=isothermal ./run_quarter_annulus_hyper_dual.sh`
- Outcome:
  - Max iteration reached: `3600`
  - NaN index (if any): `(-2,-2)`
  - First visible oscillation location: not re-measured in this run
  - Notes:
    - Replacing symmetry-side BC with extrapolation did not shift or remove the instability.
    - After this test, `theta=pi` BC was restored to `SymmetryBC`.
