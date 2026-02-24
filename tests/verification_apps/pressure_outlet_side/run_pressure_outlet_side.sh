#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OPS_DIR_DEFAULT="${ROOT_DIR}/../OPS"
PYTHON_BIN="${PYTHON_BIN:-python3}"
OPS_APP_BIN_DIR="${OPS_APP_BIN_DIR:-${OPS_DIR_DEFAULT}/ops/c/lib/gnu}"
HDF5_LIB_DIR="${HDF5_LIB_DIR:-}"

NITER_VALUE="${NITER_VALUE:-200}"
DT_VALUE="${DT_VALUE:-2.0e-4}"
NP0_VALUE="${NP0_VALUE:-81}"
NP1_VALUE="${NP1_VALUE:-41}"

cd "${CASE_DIR}"
export PYTHONPATH="${ROOT_DIR}:${PYTHONPATH:-}"

if [[ ! -d "${OPS_APP_BIN_DIR}" ]]; then
  echo "OPS library directory not found: ${OPS_APP_BIN_DIR}"
  exit 1
fi

if [[ -z "${HDF5_LIB_DIR}" ]]; then
  if [[ -d "/opt/homebrew/lib" ]]; then
    HDF5_LIB_DIR="/opt/homebrew/lib"
  elif [[ -d "/usr/lib/x86_64-linux-gnu" ]]; then
    HDF5_LIB_DIR="/usr/lib/x86_64-linux-gnu"
  else
    HDF5_LIB_DIR="/usr/lib"
  fi
fi

run_case () {
  local side="$1"
  echo "[case side=${side}] generate"
  OUTLET_SIDE="${side}" \
  NITER_VALUE="${NITER_VALUE}" \
  DT_VALUE="${DT_VALUE}" \
  NP0_VALUE="${NP0_VALUE}" \
  NP1_VALUE="${NP1_VALUE}" \
  "${PYTHON_BIN}" verify_pressure_outlet_side.py

  echo "[case side=${side}] translate"
  "${PYTHON_BIN}" "${ROOT_DIR}/../OPS/ops_translator_legacy/c/ops.py" opensbli.cpp

  echo "[case side=${side}] build"
  clang++ -O3 -fPIC -Wall -g -std=c++11 -Dgnu \
    -I"${ROOT_DIR}/../OPS/ops/c/include" \
    opensbli_ops.cpp -I. ./mpi_openmp/mpi_openmp_kernels.cpp \
    -L"${OPS_APP_BIN_DIR}" -lops_seq -lops_hdf5_seq \
    -L"${HDF5_LIB_DIR}" -lhdf5_hl -lhdf5 -lz \
    -o opensbli_seq

  echo "[case side=${side}] run"
  ./opensbli_seq | tee "run_side${side}.log"
}

# Run reference (original implementation path) and side0 path.
run_case 1
run_case 0

echo "Done. Compare run_side1.log vs run_side0.log and monitor files."
