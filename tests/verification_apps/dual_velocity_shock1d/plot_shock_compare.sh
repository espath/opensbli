#!/usr/bin/env bash
set -euo pipefail

THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${THIS_DIR}"

MA_VALUE="${MA_VALUE:-1.55}"
ML_VALUE="${ML_VALUE:-0.4}"
NP0_VALUE="${NP0_VALUE:-401}"
ma_tag="$(printf '%g' "${MA_VALUE}" | sed 's/\./p/g')"
ml_tag="$(printf '%g' "${ML_VALUE}" | sed 's/\./p/g')"
dual_file="outputs/opensbli_output_ma${ma_tag}_ml${ml_tag}_np${NP0_VALUE}.h5"
nsf_file="outputs/opensbli_output_ma${ma_tag}_ml0_np${NP0_VALUE}.h5"

# If explicit NP0 file is missing, auto-pick a mesh-tagged file if available.
if [[ ! -f "${dual_file}" ]]; then
  cand_dual="$(ls -1 outputs/opensbli_output_ma"${ma_tag}"_ml"${ml_tag}"_np*.h5 2>/dev/null | head -n 1 || true)"
  if [[ -n "${cand_dual}" ]]; then
    dual_file="${cand_dual}"
  else
    # Backward-compatible fallback for pre-np-tag output files.
    dual_file="outputs/opensbli_output_ma${ma_tag}_ml${ml_tag}.h5"
  fi
fi
if [[ ! -f "${nsf_file}" ]]; then
  cand_nsf="$(ls -1 outputs/opensbli_output_ma"${ma_tag}"_ml0_np*.h5 2>/dev/null | head -n 1 || true)"
  if [[ -n "${cand_nsf}" ]]; then
    nsf_file="${cand_nsf}"
  else
    # Backward-compatible fallback for pre-np-tag output files.
    nsf_file="outputs/opensbli_output_ma${ma_tag}_ml0.h5"
  fi
fi

if [[ ! -f "${dual_file}" ]]; then
  echo "Missing dual-velocity output: ${dual_file}"
  echo "Available for Ma=${ma_tag}:"
  ls -1 outputs/opensbli_output_ma"${ma_tag}"_ml*.h5 2>/dev/null || true
  echo "Generate it with:"
  echo "  MA_VALUE=${MA_VALUE} ML_VALUE=${ML_VALUE} ./run_shock_dual.sh"
  exit 1
fi

if [[ ! -f "${nsf_file}" ]]; then
  echo "Missing NSF output: ${nsf_file}"
  echo "Generate it with:"
  echo "  MA_VALUE=${MA_VALUE} ML_VALUE=0.0 ./run_shock_nsf.sh"
  exit 1
fi

python3 plot_compare_density_opensbli.py \
  --Ma "${MA_VALUE}" \
  --Ml "${ML_VALUE}" \
  --opensbli "${dual_file}" \
  --opensbli2 "${nsf_file}" \
  --label2 'OpenSBLI NSF ($\mathrm{M}_\ell=0$)' \
  --out "shock_compare_dual_nsf_exp_ma${ma_tag}_ml${ml_tag}.pdf"
