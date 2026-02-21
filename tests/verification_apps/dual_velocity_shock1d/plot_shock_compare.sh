#!/usr/bin/env bash
set -euo pipefail

THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${THIS_DIR}"

MA_VALUE="${MA_VALUE:-1.55}"
ML_VALUE="${ML_VALUE:-0.4}"
ma_tag="${MA_VALUE//./p}"
ml_tag="${ML_VALUE//./p}"

python3 plot_compare_density_opensbli.py \
  --Ma "${MA_VALUE}" \
  --Ml "${ML_VALUE}" \
  --opensbli "outputs/opensbli_output_ma${ma_tag}_ml${ml_tag}.h5" \
  --opensbli2 "outputs/opensbli_output_ma${ma_tag}_ml0.h5" \
  --label2 'OpenSBLI NSF ($\mathrm{M}_\ell=0$)' \
  --out "shock_compare_dual_nsf_exp_ma${ma_tag}_ml${ml_tag}.pdf"
