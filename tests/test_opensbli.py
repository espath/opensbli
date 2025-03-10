"""Integration tests to ensure the chosen test apps can be generated, translated
and compiled to target frameworks.

This module provides a set of functions for testing the generation, translation,
and compilation of test applications to target frameworks. These functions are
used to validate the functionality and compatibility of the applications with
the target frameworks.
"""

from __future__ import annotations

import glob
import logging
import os
import pathlib
import shutil
import subprocess
from argparse import ArgumentParser
from enum import Enum, auto
from functools import reduce
from typing import Dict, List


class TestModes(Enum):
    """
    Enum class for the different test modes.
    """

    ALL_TEST_CASES = auto()
    SIMPLE_TEST_CASES = auto()
    VERIFICATION_TEST_CASES = auto()


class TranslatorMode(Enum):
    """
    Enum class for the different translator modes.
    """

    LEGACY = auto()
    MODERN = auto()


# Default is to run all test cases
TEST_MODE = TestModes.ALL_TEST_CASES
SCRIPT_DIRECTORY = os.path.dirname(os.path.realpath(__file__))
VERBOSE = False
TRANSLATOR_MODE = TranslatorMode.MODERN

# The following two lists are used to define the test cases to be run. The first
# list contains the test applications, while the second list contains the
# verification applications. The verification applications are used to check the
# output of the test applications against expected values.
APP_TEST_CASES = [
    f"{SCRIPT_DIRECTORY}/../apps/wave/wave.py",
    f"{SCRIPT_DIRECTORY}/../apps/euler_wave/euler_wave.py",
    f"{SCRIPT_DIRECTORY}/../apps/shu_osher/shu_osher.py",
    f"{SCRIPT_DIRECTORY}/../apps/Sod_shock_tube/Sod_shock_tube.py",
    f"{SCRIPT_DIRECTORY}/../apps/Lax_shock_tube/Lax_shock_tube.py",
    f"{SCRIPT_DIRECTORY}/../apps/LeBlanc/LeBlanc.py",
    f"{SCRIPT_DIRECTORY}/../apps/taylor_green_vortex/taylor_green_vortex.py",
    f"{SCRIPT_DIRECTORY}/../apps/taylor_green_vortex/TGsym/TGsym.py",
    f"{SCRIPT_DIRECTORY}/../apps/viscous_shock_tube/viscous_shock_tube.py",
    f"{SCRIPT_DIRECTORY}/../apps/kelvin_helmholtz/kelvin_helmholtz.py",
    f"{SCRIPT_DIRECTORY}/../apps/inviscid_shock_reflection/inviscid_shock.py",
    f"{SCRIPT_DIRECTORY}/../apps/katzer_SBLI/katzer_SBLI.py",
    f"{SCRIPT_DIRECTORY}/../apps/channel_flow/laminar_2D/laminar_channel.py",
    f"{SCRIPT_DIRECTORY}/../apps/channel_flow/turbulent_3D/turbulent_channel.py",
    f"{SCRIPT_DIRECTORY}/../apps/channel_flow/compressible_TCF_Central/turbulent_channel.py",
    f"{SCRIPT_DIRECTORY}/../apps/channel_flow/compressible_TCF_TENO/turbulent_channel.py",
    f"{SCRIPT_DIRECTORY}/../apps/channel_flow/adiabatic_isothermal_channel/iso_adi_channel_heat_sink.py",
    f"{SCRIPT_DIRECTORY}/../apps/transitional_SBLI/transitional_SBLI.py",
    f"{SCRIPT_DIRECTORY}/../apps/cylinder/supersonic_cylinder/supersonic_cylinder.py",
    f"{SCRIPT_DIRECTORY}/../apps/compressible_taylor_green_vortex/TGV_multi_block/compressible_TGV_MB.py",
    f"{SCRIPT_DIRECTORY}/../apps/compressible_taylor_green_vortex/compressible_TGV.py",
    f"{SCRIPT_DIRECTORY}/../apps/aerofoils/multi_block/2D/airfoil_MB_2D.py",
    f"{SCRIPT_DIRECTORY}/../apps/aerofoils/multi_block/3D/transonic_MB.py",
]
VERIFICATION_TEST_CASES = [
    f"{SCRIPT_DIRECTORY}/verification_apps/euler_wave_curvilinear/verify_euler_wave_curvilinear.py",
    f"{SCRIPT_DIRECTORY}/verification_apps/sod_shock_tube/verify_sod_shock_tube.py",
    f"{SCRIPT_DIRECTORY}/verification_apps/tg_sym/verify_tg_sym.py",
    f"{SCRIPT_DIRECTORY}/verification_apps/channel_flow_laminar_2D/verify_laminar_channel.py",
    f"{SCRIPT_DIRECTORY}/verification_apps/airfoil_multiblock_2D/verify_airfoil_MB_2D.py",
]
ALL_TEST_CASES = APP_TEST_CASES + VERIFICATION_TEST_CASES


def setup_logger(log_file: str) -> logging.Logger:
    """
    Setup the file logger.

    Parameters
    ----------
    log_file : str
        The name of the log file to write to, including file path.

    Returns
    -------
    logging.Logger
        The logger object
    """
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)

    # Create a file handler and set formatting
    file_handler = logging.FileHandler(log_file, mode="w")
    file_handler.setLevel(logging.INFO)
    formatter = logging.Formatter("%(message)s")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    return logger


logger = setup_logger("./test.log")


def _log(message: str, suppress_output: bool = False) -> None:
    """
    Print a message to standard out.

    Parameters
    ----------
    message : str
        The message to be printed.
    suppress_output : bool
        Suppress output to the screen, but will be overridden by global
        verbosity setting.
    """
    logger.info(f"{message}")
    if not suppress_output or VERBOSE:
        print(message)


def _log_failure(message: str, suppress_output: bool = False) -> None:
    """
    Print a failure message in red.

    Parameters
    ----------
    message : str
        The message to be printed.
    suppress_output : bool
        Suppress output to the screen, but will be overridden by global
        verbosity setting.
    """
    logger.info(f"{message}")
    if not suppress_output or VERBOSE:
        print(f"\033[91m⨉ {message}\033[0m")


def _log_success(message: str, suppress_output: bool = False) -> None:
    """
    Print a success message in green.

    Parameters
    ----------
    message : str
        The message to be printed.
    suppress_output : bool
        Suppress output to the screen, but will be overridden by global
        verbosity setting.
    """
    logger.info(f"{message}")
    if not suppress_output or VERBOSE:
        print(f"\033[92m✓ {message}\033[0m")


def _run_process(
    commands: list, cwd: str, use_shell: bool = False
) -> Dict[int, str | None, str | None]:
    """
    Execute a subprocess with the given commands and current working directory.

    Parameters
    ----------
    commands : list of str
        The list of command-line arguments to be executed.
    cwd : str
        The path to the current working directory for the subprocess.
    use_shell : bool
        Execute the provided commands directly in the shell. Useful for when
        you need to activate a Python environment first.

    Returns
    -------
    dict
        A dictionary containing the exit code, standard output, and standard
        error.

    Notes
    -----
    This function executes the specified commands as a subprocess. If the
    `VERBOSE` global variable is set to True, it prints the standard output of
    the subprocess. If the subprocess returns a non-zero exit status, it prints
    the failed command along with the stderr output.
    """
    _log(f"- Running: `{' '.join(commands)}`", suppress_output=True)

    try:
        if use_shell:
            rc = subprocess.run(
                " ".join(commands),
                cwd=cwd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                executable="/bin/bash",
                shell=True,
            )
        else:
            rc = subprocess.run(
                commands,
                cwd=cwd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
        return {
            "return_code": rc.returncode,
            "stdout": rc.stdout.decode(),
            "stderr": rc.stderr.decode()
            if rc.stderr
            else "Process failed with no error message",
        }
    except Exception as e:
        return {"return_code": 1, "stdout": None, "stderr": str(e)}


def _activate_python_env(commands: list, env_location: str) -> list:
    """
    Prepend an activation of a Python venv to a list of commands.

    Parameters
    ----------
    commands : list
        The list of commands to prepend with the environment activation.
    env_location : str
        The file path to the Python environment to activate.

    Returns
    -------
    list
        The updated listed of commands.

    Notes
    -----
    This function exists for when you need to activate a Python environment for
    a set of commands or script. This is a fix for an issue encountered on WSL
    where CMake would not be able to activate the OPS translator venv causing
    an app build to fail. By activating the venv first, CMake finds the correct
    venv required for translation.
    """
    return ("source", f"{env_location}/bin/activate", "&&", *commands)


def _copy_cmakelists_to_app(app_dir: str) -> str:
    """
    Copy the OPS CMakeLists.txt to the app directory.

    Parameters
    ----------
    app_dir : str
        The directory containing the application to be compiled.

    Returns
    -------
    str
        The file path of the CMakeLists.txt copied to the app directory.

    Notes
    -----
    Checks for a CMakeLists.txt in the "EnvDirectory" which is the directory
    containing the OpenSBLI, OPS, Python venv and possibly HDF5 files created
    using the "CreateOpenSBLIEnv.sh" installation script. If this is not the
    setup used, then the CMakeLists.txt in the apps/ directory will be copied
    instead. This CMakeLists.txt may be out of date with the latest version of
    OPS.
    """

    path = pathlib.Path(f"{SCRIPT_DIRECTORY}/../../CMakeLists.txt")
    if path.exists():
        # need to modify this file to remove OPS mess and add a OpenSBLI build option
        content = path.read_text()
        content = "\n".join(
            [line for line in content.splitlines() if "add_subdirectory" not in line]
        )
        content += '\nBUILD_OPS_C_SAMPLE(OpenSBLI "NONE" "NONE" "NONE" "NO" "NO")\n'
        new_path = pathlib.Path(f"{app_dir}/CMakeLists.txt")
        new_path.write_text(content)
    else:
        path = pathlib.Path(f"{SCRIPT_DIRECTORY}/../apps/CMakeLists.txt")
        shutil.copyfile(str(path.absolute()), f"{app_dir}/CMakeLists.txt")

    return str(path.resolve())


def prepare_test_environment() -> List[str]:
    """
    Prepare the test environment.

    Retrieves a list of paths to test cases, as defined in the _test_cases.py
    module. Test cases are copied into an _opensbli-test-workspace directory.

    Returns
    -------
    list
        A list of paths to various test application scripts within the root
        directory.

    Notes
    -----
    The location of app scripts are assumed to be a fixed
    """
    _log("Preparing test environment")

    # Clean up and create workspace directory for running the tests
    workspace_dir = f"{SCRIPT_DIRECTORY}/_opensbli-test-workspace"
    shutil.rmtree(workspace_dir, ignore_errors=True)
    os.makedirs(workspace_dir)

    # Set the test cases depending on chosen test mode. The default is test
    # everything
    if TEST_MODE is TestModes.VERIFICATION_TEST_CASES:
        test_cases = VERIFICATION_TEST_CASES
    else:
        test_cases = ALL_TEST_CASES

    # Copy the tests defined above into the workspace directory. We are
    # retrieving the test cases from _test_cases.py
    workspace_tests = []
    for test_case in test_cases:
        old_dir = os.path.dirname(test_case)
        new_dir_base = f"{workspace_dir}/{old_dir.split(os.path.sep)[-1]}"
        _log(
            f"- Copying {os.path.basename(test_case)} to _opensbli-test-workspace/{os.path.basename(old_dir)}",
            suppress_output=True,
        )
        new_dir = new_dir_base
        counter = 2  # in case there are cases with the same directory name
        while os.path.exists(new_dir):
            new_dir = f"{new_dir_base}-{counter}"
            counter += 1

        copied = shutil.copytree(
            old_dir,
            new_dir,
        )
        workspace_tests.append(f"{copied}/{os.path.basename(test_case)}")

    return workspace_tests


def test_app_generates(app_file: str, app_dir: str) -> int:
    """
    Generate a test app with the specified file and directory.

    Parameters
    ----------
    app_file : str
        The path to the application file to be generated.
    app_dir : str
        The directory where the application will be generated.

    Returns
    -------
    int
        The return code of the process execution.

    Raises
    ------
    ValueError
        If `app_file` or `app_dir` is not provided.

    Notes
    -----
    This function generates an application using the specified file and
    directory by running a subprocess with the command "python <app_file>" in
    the specified directory. For this to work, OpenSBLI must be in PYTHONPATH.
    """
    if not app_file or not app_dir:
        raise ValueError("Both `app_file` and `app_dir` must be provided")

    rc = _run_process(("python", app_file), app_dir)
    if rc["return_code"]:
        _log_failure(f"OpenSBLI failed to generate C source code:\n{rc['stderr']}")

    return rc["return_code"]


def test_app_translates(app_dir: str) -> int:
    """
    Translate generated OpenSBLI code to OPS parallel code with the specified
    directory.

    Parameters
    ----------
    app_dir : str
        The directory where the application code is located.

    Returns
    -------
    int
        The return code of the translation process.

    Raises
    ------
    ValueError
        If `app_dir` is not provided.
    EnvironmentError
        If the OPS_TRANSLATOR environment variable is not set.

    Notes
    -----
    This function translates the generated OpenSBLI code to OPS parallel code
    using an external translator tool specified in the OPS_TRANSLATOR
    environment variable.The translator is expected to be located in the
    OPS_TRANSLATOR directory and named 'ops.py'.
    """
    if not app_dir:
        raise ValueError("app_dir must be provided.")

    translator = os.getenv("OPS_TRANSLATOR")
    if not translator:
        raise EnvironmentError("OPS_TRANSLATOR environment variable is not set.")
    commands = (
        "python",
        f"{translator}/ops.py",
        "opensbli.cpp",
    )
    if TRANSLATOR_MODE == TranslatorMode.MODERN:
        commands = _activate_python_env(
            commands, f"{os.getenv('OPS_TRANSLATOR')}/../ops_venv"
        )
    rc = _run_process(commands, app_dir, use_shell=True)
    if rc["return_code"]:
        _log_failure(
            f"Failed to translate using OPS_TRANSLATOR into parallel OPS code:\n{rc['stderr']}"
        )

    return rc["return_code"]


def test_app_cmake_build(app_dir: str) -> int:
    """
    Compile an OpenSBLI (and OPS translated) application located in the
    specified directory.

    Parameters
    ----------
    app_dir : str
        The directory containing the application to be compiled.

    Returns
    -------
    int
        The return code of the compilation process. 0 indicates success,
        non-zero values indicate errors.

    Raises
    ------
    ValueError
        If `app_dir` is not provided.

    Notes
    -----
    This function compiles an OpenSBLI application located in the specified
    directory by copying the required Makefile from the OPENSBLI_INSTALL
    directory and then executing 'make' for the parallel frameworks which
    are available.
    """
    if not app_dir:
        raise ValueError("app_dir must be provided.")

    # Copy CMake file to test directory
    # create build dir, where target binaries will be built and use CMake to
    # prepare build files and build the test app
    build_dir = f"{app_dir}/test-build"
    shutil.rmtree(build_dir, ignore_errors=True)
    os.makedirs(build_dir, exist_ok=True)

    # Prepare build files using CMake
    commands = (
        "cmake",
        "..",
        f"-DOPS_INSTALL_DIR={os.getenv('OPS_INSTALL_DIR')}",
        "-DCMAKE_BUILD_TYPE=Release",
        f"-DLEGACY_CODEGEN={'ON' if TRANSLATOR_MODE is TranslatorMode.LEGACY else 'OFF'}",
    )
    # add -DHDF5_ROOT if HDF5_INSTALL_PATH env variable is found
    if os.getenv("HDF5_INSTALL_PATH"):
        commands = (
            *commands,
            *(f"-DHDF5_ROOT={os.getenv('HDF5_INSTALL_PATH')}",),
        )
    if TRANSLATOR_MODE == TranslatorMode.MODERN:
        commands = _activate_python_env(
            commands, f"{os.getenv('OPS_TRANSLATOR')}/../ops_venv"
        )
    rc = _run_process(commands, build_dir, use_shell=True)
    if rc["return_code"]:
        shutil.rmtree(build_dir)
        _log_failure(f"CMake failed to configure properly:\n{rc['stderr']}")
        return rc["return_code"]

    # Build apps using CMake. Use only 1 core to make output clearer
    rc = _run_process(("cmake", "--build", "."), build_dir)
    if rc["return_code"]:
        _log_failure(f"Failed to translate and build apps using CMake:\n{rc['stderr']}")

    return rc["return_code"]


def test_app_output(app_dir: str) -> int:
    """
    Test the output of an application.

    Parameters
    ----------
    app_dir : str
        The directory containing the applications to be run

    Returns
    -------
    int
        Exit code indicating the overall success or failure of the tests.
        A value of 0 indicates success, while non-zero values indicate failure.

    Raises
    ------
    ValueError
        If `app_dir` is not provided.

    Notes
    -----
    This function tests the output of an application located in `app_dir`.
    It runs the application for different modes and verifies the output.

    The function first runs the application, assuming it has already been
    compiled. If the application runs successfully, it proceeds to run a
    verification script, which checks the output for similarity against expected
    values.
    """
    if not app_dir:
        raise ValueError("app_dir must be provided.")
    if not os.path.exists(f"{app_dir}/check_opensbli_output.py"):
        return 0

    build_dir = f"{app_dir}/test-build"
    globbed_apps = [
        os.path.basename(app)
        for app in sorted(glob.glob(f"{build_dir}/OpenSBLI_*[!_opencl*]"))
    ]
    _log(f"Apps found: {globbed_apps}", suppress_output=True)

    mode_rc = {}
    for app_name in globbed_apps:
        # first run the application, as it should have compiled
        shutil.move(f"{build_dir}/{app_name}", f"{app_dir}/{app_name}")
        if "mpi" in app_name:
            mode = _run_process(("mpirun", "-n", "2", f"./{app_name}"), app_dir)
        else:
            mode = _run_process((f"./{app_name}",), app_dir)
        # Don't need to keep it around, especially since it pollutes git
        os.remove(f"{app_dir}/{app_name}")
        # Skip the next step if the app fails to run completely
        if mode["return_code"]:
            mode_rc[app_name] = mode
            continue
        # then run verification script which will return non-zero if the tested
        # quantities are not similar enough
        # todo: modify to pass relative tolerance as argument  to script
        mode_rc[app_name] = _run_process(
            ("python", "check_opensbli_output.py"), app_dir
        )
    # count how many failed using reduction
    return_code = reduce(lambda x, y: x + abs(y["return_code"]), mode_rc.values(), 0)
    shutil.rmtree(build_dir, ignore_errors=True)

    # print error messages, if there are any
    if return_code:
        for name, mode in mode_rc.items():
            if mode["return_code"]:
                _log_failure(
                    f"{name} failed with error code {mode['return_code']}:\n{mode['stderr']}",
                )

    return return_code


def run_tests() -> None:
    """
    Run the testing framework.

    Raises
    ------
    EnvironmentError
        If the environment variables OPENSBLI_INSTALL or OPS_TRANSLATOR are
        not set.
    """
    # the following environment variables are all required to be set
    if not os.getenv("OPS_INSTALL_DIR"):  # todo: might be required for CMake build
        raise EnvironmentError("$OPS_INSTALL_DIR has not been set")
    if not os.getenv("OPS_TRANSLATOR"):
        raise EnvironmentError("$OPS_TRANSLATOR has not been set")

    _log(f"OPS_INSTALL_DIR : {os.getenv('OPS_INSTALL_DIR')}", suppress_output=True)
    _log(f"OPS_TRANSLATOR  : {os.getenv('OPS_TRANSLATOR')}", suppress_output=True)

    test_app_paths = prepare_test_environment()
    num_tests = len(test_app_paths)
    num_failed = 0

    for app_path in test_app_paths:
        app_name = os.path.splitext(os.path.basename(app_path))[0]
        app_file = os.path.basename(app_path)
        app_dir = os.path.dirname(app_path)
        cmake_src = _copy_cmakelists_to_app(app_dir)
        _log("-" * 80)
        _log(f"Testing: \033[1m{app_name}\033[0m")
        _log(f"Directory: {app_dir}", suppress_output=True)
        _log(f"CMakeLists.txt: {cmake_src}", suppress_output=True)

        # Test that OpenSBLI can generate the test case
        return_code = test_app_generates(app_file, app_dir)
        if return_code:
            num_failed += 1
            continue
        # Test that OPS can translate the OpenSBLI code
        return_code = test_app_translates(app_dir)
        if return_code:
            num_failed += 1
            continue
        # Test that the OPS translated code can compile
        return_code = test_app_cmake_build(app_dir)
        if return_code:
            num_failed += 1
            continue
        # Test that the compile apps run properly and produce the correct
        # output
        return_code = test_app_output(app_dir)
        if return_code:
            num_failed += 1
            continue

        _log_success("Passed all tests successfully")

    _log("-" * 80)
    if num_failed > 0:
        _log_failure(f"{num_tests - num_failed}/{num_tests} tests passed")
    else:
        _log_success(f"{num_tests - num_failed}/{num_tests} tests passed")
    _log("-" * 80)

    if num_failed == 0:
        shutil.rmtree(f"{SCRIPT_DIRECTORY}/_opensbli-test-workspace")

    return num_failed


if __name__ == "__main__":
    # Parse command line arguments
    ap = ArgumentParser()
    ap.add_argument("--verif-only", action="store_true", default=False)
    ap.add_argument("--legacy-translator", action="store_true", default=False)
    ap.add_argument("--verbose", action="store_true", default=False)
    args = ap.parse_args()

    # Set the global verbosity flag
    VERBOSE = args.verbose
    # Set the test mode based on the command line arguments
    if args.verif_only:
        TEST_MODE = TestModes.VERIFICATION_TEST_CASES
    if args.legacy_translator:
        TRANSLATOR_MODE = TranslatorMode.LEGACY

    # Run the test procedure
    _log(f"Script directory: {SCRIPT_DIRECTORY}", suppress_output=True)
    num_failed = run_tests()
    exit(num_failed)
