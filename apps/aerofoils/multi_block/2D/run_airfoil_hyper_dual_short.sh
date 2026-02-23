#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"

# shellcheck disable=SC1091
source "${ROOT_DIR}/env_ops.sh"

PYTHON_BIN="${PYTHON_BIN:-${ROOT_DIR}/espath312/bin/python3}"
export PYTHONPATH="${ROOT_DIR}:${PYTHONPATH:-}"

# Safe short-run defaults (override from shell if needed)
export NITER_VALUE="${NITER_VALUE:-2000}"
export DT_VALUE="${DT_VALUE:-5.0e-5}"
export ML_VALUE="${ML_VALUE:-0.4}"
export ENABLE_SHOCK_CAPTURING="${ENABLE_SHOCK_CAPTURING:-0}"
export ENABLE_DRP_FILTER="${ENABLE_DRP_FILTER:-0}"
export ENABLE_BINOMIAL_FILTER="${ENABLE_BINOMIAL_FILTER:-1}"
export ENABLE_NAN_CHECK="${ENABLE_NAN_CHECK:-1}"

echo "[cleanup] removing stale outputs"
(
  cd "${SCRIPT_DIR}"
  rm -f opensbli_output*.h5 metrics.h5
)

echo "[1/3] Generate OpenSBLI code (niter=${NITER_VALUE}, dt=${DT_VALUE}, Ml=${ML_VALUE})"
(
  cd "${SCRIPT_DIR}"
  "${PYTHON_BIN}" airfoil_MB_2D_hyper_dual.py
)

echo "[2/3] Build opensbli_seq"
(
  cd "${SCRIPT_DIR}"
  make clean >/dev/null 2>&1 || true
  make opensbli_seq \
    OPS_COMPILER="${OPS_COMPILER:-gnu}" \
    CXX=clang++ CC=clang MPICXX=clang++ MPICC=clang \
    OMPFLAGS="" \
    CXXFLAGS="-O3 -fPIC -Wall -g -std=c++11" \
    CCFLAGS="-O3 -std=c99 -fPIC -Wall -g"
)

echo "[3/4] Run opensbli_seq"
(
  cd "${SCRIPT_DIR}"
  ./opensbli_seq
)

echo "[4/4] Tag outputs with Ml=${ML_VALUE}"
(
  cd "${SCRIPT_DIR}"
  ML_TAG="$(printf "%s" "${ML_VALUE}" | sed 's/-/m/g; s/\./p/g')"
  shopt -s nullglob
  for f in opensbli_output*.h5; do
    cp -f "${f}" "${f%.h5}_ml${ML_TAG}.h5"
  done
  shopt -u nullglob
)
