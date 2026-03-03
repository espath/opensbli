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
NITER_VALUE="${NITER_VALUE:-200}"
NP0_VALUE="${NP0_VALUE:-}"
NP1_VALUE="${NP1_VALUE:-}"
SAVE_EVERY_VALUE="${SAVE_EVERY_VALUE:-1000}"
GRID_FILE="${GRID_FILE:-${GRID_DIR}/data_full.h5}"
PERIODIC_DIR="${PERIODIC_DIR:-0}"
VALIDATE_MESH="${VALIDATE_MESH:-1}"
SEAM_CHECK_DEPTH="${SEAM_CHECK_DEPTH:-5}"
SEAM_COORD_TOL="${SEAM_COORD_TOL:-1e-10}"
SEAM_METRIC_TOL="${SEAM_METRIC_TOL:-1e-2}"
CHECK_BC_OVERLAP="${CHECK_BC_OVERLAP:-1}"
BC_OVERLAP_MODE="${BC_OVERLAP_MODE:-report}" # report | strict

cd "${CASE_DIR}"

echo "[cleanup] removing stale outputs"
rm -f opensbli_output*.h5

if [[ ! -f "${GRID_FILE}" ]]; then
  echo "Missing mesh file: ${GRID_FILE}"
  echo "Generate it with:"
  echo "  cd ${GRID_DIR} && /Users/pmzle/Documents/dev/opensbli/espath312/bin/python3 generate_cylinder_annulus.py --nr 241 --ntheta 721 --r-inner 0.5 --r-outer 20.0 --theta0 0.0 --theta1 6.283185307179586 --beta-r 2.0 --transpose --output data_full.h5"
  exit 1
fi

ln -sf "${GRID_FILE}" data.h5

if [[ -z "${NP0_VALUE}" || -z "${NP1_VALUE}" ]]; then
  read -r MESH_NP0 MESH_NP1 < <("${PYTHON_BIN}" - <<'PY'
import h5py
with h5py.File('data.h5', 'r') as f:
    d = f['opensbliblock00']['x0_B0']
    size = d.attrs.get('size')
    if size is not None and len(size) == 2:
        print(int(size[0]), int(size[1]))
    else:
        s = d.shape
        print(s[1], s[0])
PY
)
  NP0_VALUE="${NP0_VALUE:-${MESH_NP0}}"
  NP1_VALUE="${NP1_VALUE:-${MESH_NP1}}"
fi

if [[ "${VALIDATE_MESH}" == "1" ]]; then
  echo "[precheck] Validate periodic seam/halo coordinates + metrics"
  "${PYTHON_BIN}" "${GRID_DIR}/validate_periodic_annulus_mesh.py" \
    --grid data.h5 \
    --periodic-dir "${PERIODIC_DIR}" \
    --max-depth "${SEAM_CHECK_DEPTH}" \
    --coord-tol "${SEAM_COORD_TOL}" \
    --metric-tol "${SEAM_METRIC_TOL}"
fi

echo "[1/4] Generate OpenSBLI code (Ma=${MINF_VALUE}, Ml=${ML_VALUE}, Re=${RE_VALUE}, np=(${NP0_VALUE},${NP1_VALUE}), dt=${DT_VALUE}, niter=${NITER_VALUE})"
MINF_VALUE="${MINF_VALUE}" \
ML_VALUE="${ML_VALUE}" \
RE_VALUE="${RE_VALUE}" \
DT_VALUE="${DT_VALUE}" \
NITER_VALUE="${NITER_VALUE}" \
NP0_VALUE="${NP0_VALUE}" \
NP1_VALUE="${NP1_VALUE}" \
SAVE_EVERY_VALUE="${SAVE_EVERY_VALUE}" \
PERIODIC_DIR="${PERIODIC_DIR}" \
"${PYTHON_BIN}" full_annulus_hyper_dual.py

if [[ ! -d "${OPS_APP_BIN_DIR}" ]]; then
  echo "OPS library directory not found: ${OPS_APP_BIN_DIR}"
  exit 1
fi

echo "[2/4] Translate OPS"
"${PYTHON_BIN}" "${ROOT_DIR}/../OPS/ops_translator_legacy/c/ops.py" opensbli.cpp

if [[ "${CHECK_BC_OVERLAP}" == "1" ]]; then
  echo "[precheck] Analyze BC/seam halo overlap in generated code"
  "${PYTHON_BIN}" "${CASE_DIR}/check_bc_overlap.py" \
    --ops-cpp "${CASE_DIR}/opensbli_ops.cpp" \
    --bc-exchanges "${CASE_DIR}/bc_exchanges.h" \
    --np0 "${NP0_VALUE}" \
    --np1 "${NP1_VALUE}" \
    --periodic-dir "${PERIODIC_DIR}" \
    --mode "${BC_OVERLAP_MODE}"
fi

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
