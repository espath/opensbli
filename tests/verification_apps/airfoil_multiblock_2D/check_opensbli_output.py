"""Post-processing script for checking actual value against target value for
airfoil_2d.
"""

import math
import sys

import h5py

h5_data = h5py.File("opensbli_output.h5", "r")
actual_rho_u = h5_data["opensbliblock01"]["rhou0_B1"][6, 2000]
target_rho_u = 0.6742651713484463
error_code = int(not math.isclose(actual_rho_u, target_rho_u, rel_tol=1e-6))

if error_code:
    print(f"Failed to validate actual rho_u ({actual_rho_u}) against target rho_u ({target_rho_u})", file=sys.stderr)
sys.exit(error_code)
