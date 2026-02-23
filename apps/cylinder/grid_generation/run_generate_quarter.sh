#!/usr/bin/env bash
set -euo pipefail

# Quarter annulus in x<=0, y>=0:
# theta in [pi/2, pi]
python3 generate_cylinder_annulus.py \
  --theta0 1.5707963267948966 \
  --theta1 3.141592653589793 \
  "$@"
