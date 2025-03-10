"""Post-processing script for checking actual value against target value for
a symmetrical Taylor Green vortex.

This script verifies the test output of the w velocity for cell 10,10,10
(the mid-point) is -2.2438363795512375e-08 to within a relative tolerance of
1e-9.
"""

import math
import sys

import h5py


def read_dataset(group, dataset):
    d_m = group["%s" % (dataset)].attrs['d_m']
    size = group["%s" % (dataset)].shape
    read_start = [abs(d+2) for d in d_m]
    read_end = [s-abs(d+2) for d, s in zip(d_m, size)]
    if len(read_end) == 2:
        read_data = group["%s" % (dataset)][read_start[0]:read_end[0], read_start[1]:read_end[1]]
    elif len(read_end) == 3:
        read_data = group["%s" % (dataset)][read_start[0]:read_end[0], read_start[1]:read_end[1], read_start[2]:read_end[2]]
    else:
        raise NotImplementedError("")
    return read_data


h5_data = h5py.File("opensbli_output_000001.h5", "r")
block = h5_data["opensbliblock00"]
actual_w = (read_dataset(block, "rhou2_B0") / read_dataset(block, "rho_B0"))[10, 10, 10]
target_w = -2.2438363795512375e-08

error_code = int(not math.isclose(actual_w, target_w, rel_tol=1e-9))

if error_code:
    print(f"Failed to validate actual w {actual_w} against target w ({target_w})", file=sys.stderr)
sys.exit(error_code)
