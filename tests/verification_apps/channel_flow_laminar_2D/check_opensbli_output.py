"""Post-processing script for checking actual value against target value for
channel_flow_laminar_2d.

This script verifies the test output of rho_v for cell 32, 5 is
-1.3663064964665273e-08 to within a relative tolerance of 1e-9.
"""

import math
import sys

import h5py

h5_data = h5py.File("opensbli_output.h5", "r")
actual_rho_v = h5_data["opensbliblock00"]["rhou1_B0"][32, 5]
target_rho_v = -1.3663064964665273e-08
error_code = int(not math.isclose(actual_rho_v, target_rho_v, rel_tol=1e-9))

if error_code:
    print(f"Failed to validate actual rho_v ({actual_rho_v}) against target rho_v ({target_rho_v})", file=sys.stderr)
sys.exit(error_code)
