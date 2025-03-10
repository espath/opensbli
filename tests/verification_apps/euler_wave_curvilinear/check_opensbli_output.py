"""Post-processing script for checking actual value against target value for
the euler_wave_curvilinear simulation.

This script verifies the test output of rho for cell 69,69 (the middle of the
grid) is 0.9984292346559068 to within a relative tolerance of 1e-9.
"""

import math
import sys

import h5py

h5_data = h5py.File("opensbli_output.h5", "r")
actual_rho = h5_data["opensbliblock00"]["rho_B0"][69,69]
target_rho = 0.9984292346559068
error_code = int(not math.isclose(actual_rho, target_rho, rel_tol=1e-9))

if error_code:
    print(f"Failed to validate actual rho_u ({actual_rho}) against target rho ({target_rho})", file=sys.stderr)
sys.exit(error_code)
