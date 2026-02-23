#!/usr/bin/env bash
set -euo pipefail

# Upper-half annulus in y>=0:
# theta in [0, pi]
python3 generate_cylinder_annulus.py \
  --theta0 0.0 \
  --theta1 3.141592653589793 \
  "$@"
