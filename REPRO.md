# Reproducible OpenSBLI + OPS Setup (macOS)

## 1) Bootstrap

```bash
cd /path/to/opensbli
./scripts/bootstrap_macos.sh /path/to/OPS
```

This creates:
- `espath312` virtual environment
- `env_ops.sh` with local OPS/HDF5/libclang paths
- OPS layout links expected by app CMake

## 2) Activate

```bash
source espath312/bin/activate
source env_ops.sh
```

## 3) Run verification case

```bash
./scripts/run_sod_verification.sh
```

## Notes

- Python is pinned to 3.12 for compatibility.
- `scripts/requirements-opensbli.txt` pins Python dependencies.
- Output file for Sod case is produced in:
  `tests/verification_apps/sod_shock_tube/build/opensbli_output.h5`
