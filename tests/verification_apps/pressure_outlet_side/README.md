# PressureOutlet side0 vs side1 tiny test

This is a minimal 2D Euler case to compare `PressureOutletBC` behavior on:

- `side=1` (historical/original path)
- `side=0` (new generalized path)

Setup:

- uniform Cartesian domain
- periodic boundaries in `y` to remove wall/corner complications there
- one `x` side is `Dirichlet` inlet, opposite `x` side is `PressureOutletBC`
- tiny grid and short run for quick A/B

## Run both sides

```bash
cd tests/verification_apps/pressure_outlet_side
./run_pressure_outlet_side.sh
```

Optional:

```bash
NITER_VALUE=400 DT_VALUE=1.0e-4 ./run_pressure_outlet_side.sh
```

Outputs:

- `run_side1.log`
- `run_side0.log`
- `pressure_outlet_side1_monitor.log`
- `pressure_outlet_side0_monitor.log`

If side0 handling is incorrect, you typically see earlier instability/NaN in
`run_side0.log` compared to `run_side1.log`.
