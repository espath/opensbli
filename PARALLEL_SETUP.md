# Parallel Setup (MPI) on This Machine

This document gives exact commands for **this macOS machine** to build a separate MPI-capable OPS install and use it from OpenSBLI.

Paths used below:

- OPS repo: `/Users/pmzle/Documents/dev/OPS`
- OpenSBLI repo: `/Users/pmzle/Documents/dev/opensbli`
- Python venv: `/Users/pmzle/Documents/dev/opensbli/espath312`
- MPI compilers: `/opt/homebrew/bin/mpicc`, `/opt/homebrew/bin/mpicxx`
- MPI HDF5: `/opt/homebrew/opt/hdf5-mpi`

## 1) Build a separate MPI OPS install

```bash
cd /Users/pmzle/Documents/dev/OPS
mkdir -p build-mpi install-mpi
cd build-mpi

cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=/Users/pmzle/Documents/dev/OPS/install-mpi \
  -DCMAKE_C_COMPILER=/opt/homebrew/bin/mpicc \
  -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx \
  -DBUILD_OPS_CXX=ON \
  -DBUILD_OPS_FORTRAN=OFF \
  -DBUILD_OPS_APPS=OFF \
  -DHDF5_ROOT=/opt/homebrew/opt/hdf5-mpi

cmake --build . -j8
cmake --install .
```

## 2) Add legacy translator links (needed by current OpenSBLI app flow)

```bash
mkdir -p /Users/pmzle/Documents/dev/OPS/install-mpi/translator/ops_translator_legacy
ln -sfn /Users/pmzle/Documents/dev/OPS/ops_translator_legacy/c \
  /Users/pmzle/Documents/dev/OPS/install-mpi/translator/ops_translator_legacy/c
ln -sfn /Users/pmzle/Documents/dev/OPS/ops_translator_legacy/fortran \
  /Users/pmzle/Documents/dev/OPS/install-mpi/translator/ops_translator_legacy/fortran
```

## 3) Quick sanity checks for OPS MPI install

```bash
ls -l /Users/pmzle/Documents/dev/OPS/install-mpi/lib/libops_mpi.a
ls -l /Users/pmzle/Documents/dev/OPS/install-mpi/lib/libops_hdf5_mpi.a
ls -l /Users/pmzle/Documents/dev/OPS/install-mpi/include/ops_seq.h
```

## 4) Use the MPI OPS install from OpenSBLI

```bash
cd /Users/pmzle/Documents/dev/opensbli
source /Users/pmzle/Documents/dev/opensbli/espath312/bin/activate

export OPS_INSTALL_DIR=/Users/pmzle/Documents/dev/OPS/install-mpi
export OPS_TRANSLATOR=/Users/pmzle/Documents/dev/OPS/ops_translator/ops-translator

# Keep include paths clean on macOS to avoid libc++ header-order breakage.
unset C_INCLUDE_PATH CPLUS_INCLUDE_PATH CPATH CPPFLAGS CFLAGS CXXFLAGS SDKROOT
```

## 5) Configure/build dual-velocity app (MPI toolchain)

```bash
cd /Users/pmzle/Documents/dev/opensbli/tests/verification_apps/dual_velocity_shock1d
python3 verify_dual_velocity_shock1d.py

rm -rf build
mkdir build
cd build

cmake .. \
  -DOPS_INSTALL_DIR="$OPS_INSTALL_DIR" \
  -DCMAKE_BUILD_TYPE=Release \
  -DLEGACY_CODEGEN=ON \
  -DCMAKE_C_COMPILER=/opt/homebrew/bin/mpicc \
  -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx \
  -DMPI_C_COMPILER=/opt/homebrew/bin/mpicc \
  -DMPI_CXX_COMPILER=/opt/homebrew/bin/mpicxx \
  -DMPI_CXX_SKIP_MPICXX=ON \
  -DHDF5_ROOT=/opt/homebrew/opt/hdf5-mpi
```

## 6) Build and run

Sequential target (fast validation):

```bash
cmake --build . --target OpenSBLI_seq -j2
./OpenSBLI_seq
```

MPI target:

```bash
cmake --build . --target OpenSBLI_mpi -j2
mpirun -n 4 ./OpenSBLI_mpi
```

## 7) Plot result

```bash
cd /Users/pmzle/Documents/dev/opensbli/tests/verification_apps/dual_velocity_shock1d
python3 plot_compare_density_opensbli.py --Ma 1.55 --Ml 0.4 --opensbli build/opensbli_output.h5
```

## Troubleshooting

- If you see `H5Pset_fapl_mpio` / `H5Pset_dxpl_mpio` undefined symbols:
  - You mixed non-MPI HDF5 with `libops_hdf5_mpi`.
  - Reconfigure with `-DHDF5_ROOT=/opt/homebrew/opt/hdf5-mpi` and clean `build/`.

- If you see libc++ header errors (`<cstddef> ... didn't find libc++ stddef.h`):
  - `unset C_INCLUDE_PATH CPLUS_INCLUDE_PATH CPATH CPPFLAGS CFLAGS CXXFLAGS SDKROOT`
  - then reconfigure from a clean `build/`.
