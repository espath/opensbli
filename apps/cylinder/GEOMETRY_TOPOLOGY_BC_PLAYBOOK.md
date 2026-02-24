# Cylinder Setup Playbook (Geometry, Topology, BC)

This file is the persistent reference for setting up new cylinder-like cases.
Use it before creating/patching any annulus/cylinder script.

## 1) Geometry First

Define physical geometry independently from OpenSBLI directions:

- Inner wall radius: `r_inner`
- Farfield radius: `r_outer`
- Angle range:
  - quarter: `[pi/2, pi]`
  - half: `[0, pi]`
  - full: `[0, 2*pi)` (endpoint excluded)

Mesh generator: `apps/cylinder/grid_generation/generate_cylinder_annulus.py`

## 2) Topology Choice

- `full annulus`:
  - one periodic seam in angular direction
  - sensitive to seam/corner/halo handling
- `half annulus`:
  - no periodic seam
  - two radial cut boundaries (same physical line family)
- `quarter annulus`:
  - two cut boundaries + less coverage
  - useful for smoke tests

Preferred debug order:
1. half annulus
2. quarter annulus
3. full annulus

## 3) Computational Direction Mapping (Critical)

Do not assume direction labels from memory.
Read from HDF5 metadata + array shape:

- `size = [np0, np1]` (interior points)
- dataset shape includes halos.

For generated non-transposed annulus grids, do not assume this globally.
Different writers/array layouts may swap index order.

Example verified from `apps/cylinder/grid_generation/data_quarter.h5`:
- `dir0` maps to angular cut lines (`theta` edges)
- `dir1` maps to radial arcs (`r` edges)

For transposed files, mapping may invert. Verify before assigning BC.

## 4) BC Assignment by Direction (Annulus Convention)

If `dir0 = radial`, `dir1 = angular`:

- `dir0, side0` = inner wall (isothermal wall)
- `dir0, side1` = outer farfield (Dirichlet freestream)

Then angular sides depend on domain:
- full annulus: periodic on `dir1`
- half annulus: cut lines on `dir1`
- quarter annulus: cut lines on `dir1` and one additional cut on `dir0` depending on parametrization

## 5) Cut-Line BC Strategy (Half/Quarter)

For upper-half annulus (`y>=0`), angular cuts represent symmetry planes.
Target BC: symmetry.

Current blocker:
- `SymmetryBC` currently triggers a codegen stencil-name bug in this setup.
- Temporary workaround BCs (Dirichlet/Extrapolation on cuts) are for debugging only.

## 6) Metric Spacing Consistency

`Delta0block0`, `Delta1block0` must match computational directions, not physical labels:

- radial spacing: `(r_outer-r_inner)/(nr-1)`
- angular spacing:
  - half/quarter closed span: `span/(ntheta-1)`
  - full periodic: `span/(ntheta)`

Common failure mode:
- Geometry/BC map looks correct, but `Delta0/Delta1` are assigned using the opposite
  direction convention. This distorts metric scaling and can break symmetry/outlet behavior.

## 7) Minimum Smoke Checklist

Before long runs:

1. mesh metadata sanity:
   - shape and `size` consistent
2. Jacobian sanity on interior:
   - finite, no sign/zero pathologies in interior
3. run `niter=1` with writes enabled
4. check non-finite counts in interior arrays
5. only then run `niter=10/50/100`

## 8) Non-Negotiable Rule

Do not debug physics formulation here (NS vs dual/hyper-dual) for annulus setup failures.
Treat physics as validated and debug mesh/topology/BC pipeline.
