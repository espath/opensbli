"""Post-processing script for checking actual value against target value for
sod shock tube.

This script verifies the test output of rho_u for cell 100 (the mid-point) is
0.10063031315938353 to within a relative tolerance of 1e-5.
"""

import math
import sys

import h5py

h5_data = h5py.File("opensbli_output.h5", "r")
actual_rho_u = h5_data["opensbliblock00"]["rhou0_B0"][100]
target_rho_u = 0.1020821132791397
error_code = int(not math.isclose(actual_rho_u, target_rho_u, rel_tol=1e-5))

if error_code:
    print(f"Failed to validate actual rho_u ({actual_rho_u}) against target rho_u ({target_rho_u})", file=sys.stderr)
sys.exit(error_code)
