#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GRID_DIR="${ROOT_DIR}/apps/cylinder/grid_generation"
OPS_DIR_DEFAULT="${ROOT_DIR}/../OPS"
PYTHON_BIN="${PYTHON_BIN:-python3}"
OPS_APP_BIN_DIR="${OPS_APP_BIN_DIR:-${OPS_DIR_DEFAULT}/ops/c/lib/gnu}"
HDF5_LIB_DIR="${HDF5_LIB_DIR:-}"

MINF_VALUE="${MINF_VALUE:-2.0}"
RE_VALUE="${RE_VALUE:-2000}"
DT_VALUE="${DT_VALUE:-2.0e-5}"
NP0_VALUE="${NP0_VALUE:-241}"
NP1_VALUE="${NP1_VALUE:-361}"
SAVE_EVERY_VALUE="${SAVE_EVERY_VALUE:-1000}"
ITER_LIST="${ITER_LIST:-10 50 100}"
PARITY_TOL="${PARITY_TOL:-1e-12}"
GRID_FILE="${GRID_FILE:-${GRID_DIR}/data_upper_half.h5}"

cd "${CASE_DIR}"

if [[ ! -f "${GRID_FILE}" ]]; then
  echo "Missing mesh file: ${GRID_FILE}"
  exit 1
fi
ln -sf "${GRID_FILE}" data.h5

if [[ -z "${HDF5_LIB_DIR}" ]]; then
  if [[ -d "/opt/homebrew/lib" ]]; then
    HDF5_LIB_DIR="/opt/homebrew/lib"
  elif [[ -d "/usr/lib/x86_64-linux-gnu" ]]; then
    HDF5_LIB_DIR="/usr/lib/x86_64-linux-gnu"
  else
    HDF5_LIB_DIR="/usr/lib"
  fi
fi

if [[ ! -d "${OPS_APP_BIN_DIR}" ]]; then
  echo "OPS library directory not found: ${OPS_APP_BIN_DIR}"
  exit 1
fi

run_case() {
  local case_py="$1"
  local tag="$2"
  local niter="$3"
  local ml="$4"

  rm -f opensbli_output*.h5 opensbli_ops.cpp OpenSBLI_seq

  MINF_VALUE="${MINF_VALUE}" \
  RE_VALUE="${RE_VALUE}" \
  DT_VALUE="${DT_VALUE}" \
  NITER_VALUE="${niter}" \
  NP0_VALUE="${NP0_VALUE}" \
  NP1_VALUE="${NP1_VALUE}" \
  SAVE_EVERY_VALUE="${SAVE_EVERY_VALUE}" \
  ML_VALUE="${ml}" \
  ENABLE_WENO=0 ENABLE_TVD=0 ENABLE_NAN_CHECK=0 \
  "${PYTHON_BIN}" "${case_py}"

  "${PYTHON_BIN}" "${ROOT_DIR}/../OPS/ops_translator_legacy/c/ops.py" opensbli.cpp

  clang++ -O3 -fPIC -Wall -g -std=c++11 -Dgnu \
    -I"${ROOT_DIR}/../OPS/ops/c/include" \
    opensbli_ops.cpp -I. ./mpi_openmp/mpi_openmp_kernels.cpp \
    -L"${OPS_APP_BIN_DIR}" -lops_seq -lops_hdf5_seq \
    -L"${HDF5_LIB_DIR}" -lhdf5_hl -lhdf5 -lz \
    -o OpenSBLI_seq

  ./OpenSBLI_seq > "run_${tag}_n${niter}.log" 2>&1

  cp -f opensbli_output.h5 "${tag}_n${niter}.h5"
  if [[ -f opensbli_output_000001.h5 ]]; then
    cp -f opensbli_output_000001.h5 "${tag}_n${niter}_000001.h5"
  fi
}

compare_case() {
  local niter="$1"
  "${PYTHON_BIN}" - <<PY
import h5py, numpy as np, sys
tol = float("${PARITY_TOL}")
f_ns = h5py.File("ns_ml0_n${niter}.h5", "r")
f_hd = h5py.File("hyper_ml0_n${niter}.h5", "r")

def interior(arr):
    dm = arr.attrs.get("d_m")
    sz = arr.attrs.get("size")
    if dm is None or sz is None:
        return np.array(arr)
    i0 = -int(dm[0]); j0 = -int(dm[1])
    ni = int(sz[0]); nj = int(sz[1])
    data = np.array(arr)
    # Common OpenSBLI layout in these files is [np1+halos, np0+halos]
    return data[j0:j0+nj, i0:i0+ni]

fields = ["rho_B0", "rhou0_B0", "rhou1_B0", "rhoE_B0"]
maxdiff = 0.0
nonfinite_total = 0
for fld in fields:
    a = interior(f_ns["opensbliblock00"][fld])
    b = interior(f_hd["opensbliblock00"][fld])
    fa = np.isfinite(a)
    fb = np.isfinite(b)
    ff = fa & fb
    nf = int(np.count_nonzero(~ff))
    nonfinite_total += nf
    if np.any(ff):
        d = float(np.max(np.abs(a[ff] - b[ff])))
        maxdiff = max(maxdiff, d)
        print(f"{fld}: max|diff|={d:.3e}, nonfinite={nf}")
    else:
        print(f"{fld}: max|diff|=nan, nonfinite={nf}")
print(f"overall max|diff|={maxdiff:.3e}")
print(f"overall nonfinite={nonfinite_total}")
if nonfinite_total > 0 or (not np.isfinite(maxdiff)) or maxdiff > tol:
    sys.exit(2)
PY
}

echo "Half-annulus NS vs hyper-dual(Ml=0) parity matrix"
echo "ITER_LIST=${ITER_LIST} RE=${RE_VALUE} Ma=${MINF_VALUE} dt=${DT_VALUE} tol=${PARITY_TOL}"

for niter in ${ITER_LIST}; do
  echo
  echo "=== niter=${niter} ==="
  echo "[1/3] NS baseline"
  run_case half_annulus_ns_baseline.py ns_ml0 "${niter}" "0.0"
  echo "[2/3] Hyper-dual Ml=0"
  run_case half_annulus_hyper_dual.py hyper_ml0 "${niter}" "0.0"
  echo "[3/3] Compare interior fields"
  compare_case "${niter}"
done

echo
echo "Parity checks passed for all requested iterations."
