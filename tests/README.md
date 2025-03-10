# Testing OpenSBLI

This directory contains integration tests, designed to verify that OpenSBLI and
its dependencies have been installed correctly and to help prevent regression
during development. There are two sets of test cases.

In the first set of tests, each app in the OpenSBLI app directory is generated,
translated and compiled. These tests will cover the validity of OpenSBLI's code
generation for a wide range of geometries, solvers and other features. It should
also confirm if OPS and OpenSBLI are correctly installed and working together in
harmony. These tests come from the `apps/` directory in the root of the
OpenSBLI directory.

In the second sets of tests, a smaller set of cases are generated, translated,
compiled and then run with their output checked against verification data. These
tests are designed to created an "end-to-end" test scenario, to ensure the
entire pipeline works as intended. Crucially, these sets of tests will check the
run time behaviour and output of each test case, to help prevent regression
creeping in from changes to OpenSBLI, or when changing the version of
dependencies such as when updating OPS. These tests come from the
`verification_apps/` directory in this directory.

## Pre-requisites

You need OPS and OpenSBLI installed, and the following environment
variables:

```bash
PYTHONPATH=$PYTHONPATH:/path/to/opensbli
OPS_INSTALL_DIR=/path/to/OPS-INSTALL                  # contains: lib/ bin/ include/
OPS_TRANSLATOR=$OPS_INSTALL_DIR/bin/ops_translator/c  # contains: e.g. ops.py
HDF5_INSTALL_PATH=/usr/local/hdf5                     # optional, depending on where HDF5 is installed
```

OpenSBLI apps are built using CMake, therefore you will need your environment
set up to use a version of OPS built using the CMake installation method.

The optional `HDF5_INSTALL_PATH` should be used if CMake cannot find the HDF5
library and include files. This usually happens if you have built HDF5 from
source and installed them to a directory not in the regular library paths. If
HDF5 cannot by find, this usually turns up in the build/compilation stage
instead of the configuration stage.

## Running tests

Tests can be run using the python script `test_opensbli.py`, e.g.,

```bash
$ pwd
/home/.../OpenSBLI/OpenSBLI/tests
$ ls
apps/ README.md test_opensbli.py
$ python test_opensbli.py
--------------------------------------------------------------------------------
Testing: test_laminar_channel
✓ Passed all tests successfully
--------------------------------------------------------------------------------
Testing: test_euler_wave_curvilinear
✓ Passed all tests successfully
--------------------------------------------------------------------------------
Testing: test_sod_shock_tube
⨉ OpenSBLI failed with error code 1:
    Failed to validate actual rho_u (-3.427986400994968e-07) against target
    rho_u (-3.5936425368534085e-07)
--------------------------------------------------------------------------------
Testing: test_tg_sym
✓ Passed all tests successfully
--------------------------------------------------------------------------------
⨉ 3/4 tests passed
--------------------------------------------------------------------------------
```

Additional options are available and can be found using the `--help` flag. More
verbose output is written to `test.log`:

```bash
$ cat test.log
Script directory: /Users/saultyevil/srsg-projects/current-projects/opensbli/opensbli_jaxa/tests
--------------------------------------------------------------------------------
Testing: test_laminar_channel
/Users/saultyevil/srsg-projects/current-projects/opensbli/opensbli_jaxa/tests/apps/channel_flow_laminar_2D
- Running `python test_laminar_channel.py`
- Running `python $OPS_TRANSLATOR/ops.py opensbli.cpp`
- Running `cmake .. -DOPS_INSTALL_DIR=$OPS_INSTALL_DIR -DCMAKE_BUILD_TYPE=Release -DHDF5_ROOT=$HDF5_INSTALL_PATH`
- Running `cmake --build .`
Apps found: ['OpenSBLI_mpi', 'OpenSBLI_mpi_dev', 'OpenSBLI_seq', 'OpenSBLI_seq_dev']
- Running `mpirun -n 2 ./OpenSBLI_mpi`
- Running `python check_opensbli_output.py`
- Running `mpirun -n 2 ./OpenSBLI_mpi_dev`
- Running `python check_opensbli_output.py`
- Running `./OpenSBLI_seq`
- Running `python check_opensbli_output.py`
- Running `./OpenSBLI_seq_dev`
- Running `python check_opensbli_output.py`
Passed all tests successfully
--------------------------------------------------------------------------------
...
```

## Adding more tests

To add more test cases, you will need to update the relevent list of test cases
at the top of `test_opensbli.py`.

To add more "simple" test cases, add the absolute path to the OpenSBLI script
(e.g. sod_shock_tube.py) to the `APP_TEST_CASES` list.

To add more verification tests, create a directory containing two files:
`verify_*.py` and `check_opensbli_output.py`. The file `test_*.py` should be an
OpenSBLI script which generates an app. The script `check_opensbli_output.py`
should be a script which checks a value from the output files and compares it to
verification data. This script should return 0 if the app output is the
same/similar, and a non-zero code when the comparison fails. Then add the
absolute path to the OpenSBLI script to the `VERIFICATION_TEST_CASES` list.
