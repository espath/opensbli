# Cylinder Grid Generation

This folder now has two functional paths:

- `read_grid.py`: imports legacy `CYL2D.dat` and writes OpenSBLI `data.h5`.
- `generate_cylinder_annulus.py`: generates a structured annulus grid directly and writes OpenSBLI `data.h5`.

## Domain presets

Quarter annulus (`x<=0, y>=0`, i.e. `theta in [pi/2, pi]`):

```bash
cd apps/cylinder/grid_generation
./run_generate_quarter.sh --nr 241 --ntheta 241 --r-inner 0.5 --r-outer 20.0 --output data_quarter.h5
```

Upper-half annulus (`y>=0`, i.e. `theta in [0, pi]`):

```bash
cd apps/cylinder/grid_generation
./run_generate_upper_half.sh --nr 241 --ntheta 361 --r-inner 0.5 --r-outer 20.0 --output data_upper_half.h5
```

## Custom options

You can call the base generator directly:

```bash
python3 generate_cylinder_annulus.py --help
```

Useful parameters:

- `--nr`, `--ntheta`: interior resolution
- `--r-inner`, `--r-outer`: inner and outer radius
- `--beta-r`: radial clustering toward cylinder wall
- `--nhalo`: halo width (default `5`)

## About the third (blunt-body benchmark) geometry

For that case, we should first lock the exact body profile and target block topology (O-grid/C-grid/multi-block). Once you provide the benchmark geometry definition, we can add a matching generator in this same folder.
