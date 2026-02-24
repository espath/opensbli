# Quarter Annulus Hyper-Dual Cylinder Case

Single-block hyper dual-velocity setup on the quarter-annulus mesh generated in `apps/cylinder/grid_generation`.

## Defaults

- `Ma = 2.0` (`MINF_VALUE`)
- `Ml = 0.0` (`ML_VALUE`)
- `Re = 1.0e4` (`RE_VALUE`)
- `dt = 1.0e-5` (`DT_VALUE`)
- `niter = 50000` (`NITER_VALUE`)

`Re=1e4` is intentionally lower than the airfoil default (`5e4`) because this mesh is much coarser.

## Boundary mapping (quarter mesh)

For `theta in [pi/2, pi]`:

- `dir0, side0`: `theta=pi/2` radial edge (`x=0, y>0`) -> pressure outlet
- `dir0, side1`: `theta=pi` radial edge (`x<0, y=0`) -> slip/symmetry plane (`InviscidWallBC`)
- `dir1, side0`: inner arc (`r=r_inner`) -> isothermal no-slip wall
- `dir1, side1`: outer arc (`r=r_outer`) -> farfield Dirichlet

## Direction/Metric mapping (critical)

For this quarter mesh file layout (`data.h5`), indexing directions are:

- `dir0` = angular (`theta`)
- `dir1` = radial (`r`)

So spacing constants must be:

- `Delta0block0 = (pi/2)/(block0np0-1)` (angular step)
- `Delta1block0 = (r_outer-r_inner)/(block0np1-1)` (radial step)

Using the opposite assignment corrupts metric scaling and can break BC behavior.

## Symmetry-plane implementation detail

This case uses `InviscidWallBC` on `theta=pi` instead of `SymmetryBC`.
Reason: in OpenSBLI, `SymmetryBC` reflects halo values, while `InviscidWallBC`
also enforces the boundary plane state (zero normal velocity) explicitly.
For this quarter-annulus setup, that is more robust for keeping `rhou1` near zero
at the symmetry edge.

## Run

From repo root:

```bash
cd apps/cylinder/quarter_annulus_hyper_dual
./run_quarter_annulus_hyper_dual.sh
```

Override parameters as needed, e.g.:

```bash
MINF_VALUE=2.0 ML_VALUE=0.0 RE_VALUE=8000 NITER_VALUE=20000 ./run_quarter_annulus_hyper_dual.sh
```

## Note on half-annulus usefulness

Half-annulus (`theta in [0, pi]`) is useful for left-to-right flow because left and right are separate edges:

- left inflow: `theta=pi`
- right outflow: `theta=0`

You do **not** need to place inflow and outflow on the same edge in that topology.
