""" Script to generate all of the current OpenSBLI test cases."""
import os, subprocess, shutil
# List of the current applications
cases = {
'/wave/'    : 'wave.py',
'/cylinder/incompressible_cylinder/' : 'cylinder_central.py',
'/euler_wave/'  : 'euler_wave.py',
'/euler_wave_curvilinear/' : 'euler_wave_curvilinear.py',
'/shu_osher/'   : 'shu_osher.py',
'/Sod_shock_tube/'  : 'Sod_shock_tube.py',
'/Lax_shock_tube/'  : 'Lax_shock_tube.py',
'/LeBlanc/' : 'LeBlanc.py',
'/taylor_green_vortex/' : 'taylor_green_vortex.py',
'/taylor_green_vortex/TGsym/'   : 'TG_IsoT.py',
'/taylor_green_vortex/TGsym/'   : 'TGsym.py',
'/taylor_green_vortex/inviscid/': 'inviscid_taylor_green_vortex.py',
'/viscous_shock_tube/'  : 'viscous_shock_tube.py',
'/kelvin_helmholtz/'    : 'kelvin_helmholtz.py',
'/inviscid_shock_reflection/'   : 'inviscid_shock.py',
'/katzer_SBLI/' : 'katzer_SBLI.py',
'/channel_flow/laminar_2D/' : 'laminar_channel.py',
'/channel_flow/compressible_TCF_Central/'   : 'turbulent_channel.py',
'/channel_flow/compressible_TCF_TENO/'  : 'turbulent_channel.py',
'/channel_flow/adiabatic_isothermal_channel'    : 'iso_adi_channel_heat_sink.py',
'/transitional_SBLI/'   : 'transitional_SBLI.py',
'/compressible_taylor_green_vortex/'    : 'compressible_TGV.py',
'/aerofoils/multi_block/2D/'    : 'airfoil_MB_2D.py',
'/aerofoils/multi_block/3D/'    : 'transonic_MB.py',
'/turbulent_counter_flow/'      : 'turbulent_counter_flow.py',
'/vortex_core/'     : 'vortex_core.py',
'/mixing_layer/'   : 'mixlay_central_scalar.py',
'/flat_plate_transition/' : 'flat_plate_transition.py',
}


directories = [x for x in cases.keys()]
file_names = [cases[x] for x in directories]

assert len(directories) == len(file_names)
print('\33[4m' + "Found %d OpenSBLI applications." % len(file_names) + '\033[0m')
# Current working directory
owd = os.getcwd()
# Optional diff between the generated codes
check_diff = False
generate = True
compile_test = True
OPS_translator_path = '~/software/OPS/ops_translator/c/ops.py'
if check_diff:
    import difflib
    # Set a directory containing previously generated C codes
    old_code_dir = os.environ['two'] + 'apps/'

# Set the version of Python to call for testing
python_command = 'python3'
#python_command = 'python3.8'

with open(os.devnull, 'w') as devnull:

    for fname, directory in zip(file_names, directories):
        print("Generating the %s application." % (directory+fname))
        output_code = subprocess.call(["{} {}".format(python_command, fname)], shell=True, cwd=owd+directory, stdout=devnull)
        if output_code == 0:
            print('\33[92m' + "%s generated successfully." % fname + '\033[0m')
            # Compare the output code to a previously generated one
            if check_diff:
                file1, file2 = old_code_dir + directory + 'opensbli.cpp', owd + directory + 'opensbli.cpp'
                text1, text2 = open(file1).readlines(), open(file2).readlines()
                for line in difflib.unified_diff(text1, text2):
                    print(line)
            if generate:
                output_code = subprocess.call(["{} {} opensbli.cpp".format(python_command, OPS_translator_path)], shell=True, cwd=owd+directory, stdout=devnull)
                if output_code == 0:
                    print('\33[92m' + "%s translated successfully." % fname + '\033[0m')
            if compile_test:
                if not os.path.isfile('Makefile'):
                    print('\33[91m' + "No Makefile for case: {}".format(directory) + '\033[0m')
                try:
                    proc = subprocess.Popen(["make -B opensbli_openmp"], shell=True, cwd=owd+directory, stdout=devnull, stderr=subprocess.PIPE)
                    for line in proc.stderr:
                        if "error" in str(line) and "linker command failed" not in str(line):
                            print('\33[91m' + "Compilation error: {}".format(line) + '\033[0m')
                except:
                        print("Compile failed.")
        else:
            print('\33[91m' + "Generation of %s%s has failed." % (fname, directory) + '\033[0m')
