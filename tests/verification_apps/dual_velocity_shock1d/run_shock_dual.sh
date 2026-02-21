#!/usr/bin/env bash
set -euo pipefail

THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${THIS_DIR}/../../.." && pwd)"

cd "${REPO_ROOT}"
source espath312/bin/activate
source env_ops.sh
cd "${THIS_DIR}"

MA_VALUE="${MA_VALUE:-9.0}"
ML_VALUE="${ML_VALUE:-1.0}"
TVD_KAPPA_VALUE="${TVD_KAPPA_VALUE:-1.0}"
TVD_DELTA_VALUE="${TVD_DELTA_VALUE:-0.5}"
TVD_EPS_VALUE="${TVD_EPS_VALUE:-1e-8}"
BUILD_DIR="${BUILD_DIR:-build}"

echo "[1/4] Generate OpenSBLI code (Ma=${MA_VALUE}, Ml=${ML_VALUE})"
MA_VALUE="${MA_VALUE}" \
ML_VALUE="${ML_VALUE}" \
TVD_KAPPA_VALUE="${TVD_KAPPA_VALUE}" \
TVD_DELTA_VALUE="${TVD_DELTA_VALUE}" \
TVD_EPS_VALUE="${TVD_EPS_VALUE}" \
python3 verify_dual_velocity_shock1d.py

echo "[2/4] Configure CMake in ${BUILD_DIR}"
mkdir -p "${BUILD_DIR}"
# Hard-clean only configure/codegen artifacts inside BUILD_DIR.
rm -rf "${BUILD_DIR}/CMakeFiles" \
       "${BUILD_DIR}/CMakeCache.txt" \
       "${BUILD_DIR}/Makefile" \
       "${BUILD_DIR}/cmake_install.cmake" \
       "${BUILD_DIR}/tmp"
cmake -S . -B "${BUILD_DIR}" -DOPS_INSTALL_DIR="${OPS_INSTALL_DIR}" -DCMAKE_BUILD_TYPE=Release -DLEGACY_CODEGEN=ON

echo "[3/4] Build OpenSBLI_seq"
cmake --build "${BUILD_DIR}" --target OpenSBLI_seq -j2

echo "[4/4] Run solver"
(
  cd "${BUILD_DIR}"
  ./OpenSBLI_seq
)
