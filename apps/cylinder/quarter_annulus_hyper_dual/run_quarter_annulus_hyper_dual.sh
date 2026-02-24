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
ML_VALUE="${ML_VALUE:-0.0}"
RE_VALUE="${RE_VALUE:-1.0e4}"
DT_VALUE="${DT_VALUE:-1.0e-5}"
NITER_VALUE="${NITER_VALUE:-50000}"
NP0_VALUE="${NP0_VALUE:-241}"
NP1_VALUE="${NP1_VALUE:-241}"
SAVE_EVERY_VALUE="${SAVE_EVERY_VALUE:-1000}"

cd "${CASE_DIR}"

if [[ ! -f "${GRID_DIR}/data_quarter.h5" ]]; then
  echo "Missing quarter mesh: ${GRID_DIR}/data_quarter.h5"
  echo "Generate it first in apps/cylinder/grid_generation"
  exit 1
fi

ln -sf "${GRID_DIR}/data_quarter.h5" data.h5

echo "[1/4] Generate OpenSBLI code (Ma=${MINF_VALUE}, Ml=${ML_VALUE}, Re=${RE_VALUE}, np=(${NP0_VALUE},${NP1_VALUE}), dt=${DT_VALUE}, niter=${NITER_VALUE})"
MINF_VALUE="${MINF_VALUE}" \
ML_VALUE="${ML_VALUE}" \
RE_VALUE="${RE_VALUE}" \
DT_VALUE="${DT_VALUE}" \
NITER_VALUE="${NITER_VALUE}" \
NP0_VALUE="${NP0_VALUE}" \
NP1_VALUE="${NP1_VALUE}" \
SAVE_EVERY_VALUE="${SAVE_EVERY_VALUE}" \
"${PYTHON_BIN}" quarter_annulus_hyper_dual.py

if [[ ! -d "${OPS_APP_BIN_DIR}" ]]; then
  echo "OPS library directory not found: ${OPS_APP_BIN_DIR}"
  exit 1
fi

echo "[2/4] Translate OPS"
"${PYTHON_BIN}" "${ROOT_DIR}/../OPS/ops_translator_legacy/c/ops.py" opensbli.cpp

echo "[3/4] Build OpenSBLI_seq"
if [[ -z "${HDF5_LIB_DIR}" ]]; then
  if [[ -d "/opt/homebrew/lib" ]]; then
    HDF5_LIB_DIR="/opt/homebrew/lib"
  elif [[ -d "/usr/lib/x86_64-linux-gnu" ]]; then
    HDF5_LIB_DIR="/usr/lib/x86_64-linux-gnu"
  else
    HDF5_LIB_DIR="/usr/lib"
  fi
fi

clang++ -O3 -fPIC -Wall -g -std=c++11 -Dgnu \
  -I"${ROOT_DIR}/../OPS/ops/c/include" \
  opensbli_ops.cpp -I. ./mpi_openmp/mpi_openmp_kernels.cpp \
  -L"${OPS_APP_BIN_DIR}" -lops_seq -lops_hdf5_seq \
  -L"${HDF5_LIB_DIR}" -lhdf5_hl -lhdf5 -lz \
  -o OpenSBLI_seq

echo "[4/4] Run solver"
./OpenSBLI_seq
