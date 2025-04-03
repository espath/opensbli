# How to use workflow scripts

> These set scripts are intended to automatically set up an environment both using and developing OpenSBLI. However, they could also be utilised separately for installing [OPS](https://github.com/OP-DSL/OPS) library, HDF5 library, and create an independent Python environment.

- Tools needed for scripts.
  These workflow scripts will require the following tools pre-installed on a Linux(Fedora/Ubuntu)/WSL machine:

  - unzip
  - wget
  - sed
  - CMake (**version>=3.18**)
  - Python (**version>=3.8**)

- Follow the next steps to check which of the above-mentioned packages are already installed on your **existing** Linux machine or which ones came along with a **fresh** installation of a Linux/WSL machine.

   - Please check which one of the above packages is already installed. This could be checked using the "which" command. E.g., **"which unzip"**, **"which wget"**, **"which sed"**. If some packages are missing, please install them using, e.g., **sudo apt install unzip**. Also, if you are using **Fedora**, you might have a different package manager such as **dnf**; in that case, use **sudo dnf install unzip**.

   - Also, check and ensure that the versions of the cmake and python packages meet the above-mentioned requirements, using the "--version" option. e.g. **"cmake --version"**, **"python3 --version"**, etc.

   - Please check if **Python >= 3.8** installation is **present** or not using the command, **python3 --version**. This should point to a Python version which is **>=3.8**. You could also check using the command "**which python3**". This should point to **/usr/bin/python3**. You can further ensure that python3 really points to python **>=3.8** using the list command **"ll /usr/bin/python3"** and should show a soft link pointing to Python >= 3.8, e.g., if you have Python 3.12 installed then you should see, **/usr/bin/python3 -> python3.12**.

   - Also, please check if the corresponding python **venv** module is installed or not, e.g., using "**python3 -m venv --help**". If this returns a message like: "/usr/bin/python3: No module named venv". This means that this package is required to be installed separately using, e.g., on a Ubuntu machine, "**sudo apt install python3.12-venv**". You may need to use sudo if you are not the root user. **venv** module is used to create isolated virtual environments for Python projects. Without this package, commands like "python3.12 -m venv myenv" may not work.

   - If **Python >= 3.8** installation is **missing**, you could install it using the apt install command, e.g., **apt install python3.12 python3.12-venv**. If you are not a root user, please use **sudo**, e.g., **sudo apt install python3.12 python3.12-venv**. Notice that this will install both "Python 3.12" and the corresponding 3.12 virtual environment module (python3.12-venv) to help create virtual environments at the same time, as will be required later. Also, if you are using Fedora, you might have a different package manager such as **dnf**; in that case, use **sudo dnf install python3.12 python3.12-venv**.

   - On HPC machines/clusters, e.g. **Archer2**, we load some of the appropriate modules using the module load command. Please check the script named "**_InstallOPS.sh_**". This script contains already tested appropriate modules which work on Archer2. If you want to use these scripts on any new HPC system, you would need to first find the appropriate modules (for python3, cmake, hdf5, GCC (GNU) compiler for cpp, etc., see **InstallOPS.sh**) and then add them to the **InstallOPS.sh** script manually.


- Once the above packages are installed, we can create the **OpenSBLI/OPS software environment**, where we can generate, translate and compile the codes

  - Step 1 Create the environment

    Copy CreateOpenSBLIEnv.sh into a directory, make it executable and run, for example, to create an environment including Python3/2, OPS.

    ```bash
    wget -c https://raw.githubusercontent.com/opensbli/environment_and_run_scripts/refs/heads/main/CreateOpenSBLIEnv.sh
    chmod a+x CreateOpenSBLIEnv.sh
    ```

    Now, the **CreateOpenSBLIEnv.sh** script can be used to install the full **OpenSBLI/OPS software** environment. A typical command looks as follows:
    
    ```bash
    ./CreateOpenSBLIEnv.sh -d /home/username/tmp/OpenSBLIEnv
    ```    
    
    The first argument of the script above, specified through a "**-d**" flag, is for the directory where you want to install the full environment. Here, the complete installation of the **OpenSBLI/OPS software** environment will be done within the directory called **OpenSBLIEnv** inside /home/username/tmp/. The path to this directory ****must**** use an **absolute path** as shown in the command,i.e., **/home/username/tmp/OpenSBLIEnv**. Also, you must ensure that there is no directory with the same name preexisting at this path, e.g., from previous failed attempts. Otherwise, the installation will simply exit without performing any installation. Please either delete the existing faulty installation (rm -rf /home/username/tmp/OpenSBLIEnv) or give it a different name, e.g., **OpenSBLIEnv_2**.

    By default, the OPS library, the Python environment, and the OpenSBLI framework are installed inside the directory **/home/username/tmp/OpenSBLIEnv**. If required, a local HDF5 can also be installed.
    
    During the creation, the script will assume all dependencies are ready by default. But we also provide help on installing these dependencies for a few Linux distributions and HPC clusters. This can be done using an additional flag. i.e., "-m". 

    ```bash
    ./CreateOpenSBLIEnv.sh -d /home/username/tmp/OpenSBLIEnv -m Ubuntu
    ```
    So far, the script-based installation has been tested on: Linux (Ubuntu/Fedora), ARCHER2, IRIDIS5, CIRRUS (GPU machine, going out of service in April 2025). One can use the machine flag **-m** with the machine name, as shown in the Ubuntu example above. Machine names such as **Ubuntu, Fedora, ARCHER2, and IRIDIS5** should work out of the box (names are **case-sensitive**; use as mentioned).
    
    Alternatively, one could replace /home/username/ with **~/**, which is a shortcut for **/home/username/**, as below
    ```bash
    ./CreateOpenSBLIEnv.sh -d ~/tmp/OpenSBLIEnv -m Ubuntu
    ```
  - Step 2 Generate C/C++ codes for an application

    By default, Generate.sh will set the directory where it is to be the environment directory. It will activate the Python environment in the environment and translate the specified OpenSBLI Python source code.

    ```bash
    # assuming under the environment directory
    cd OpenSBLI/apps/transitional_SBLI/
    ../../../Generate.sh transitional_SBLI.py
    ```
  - Step 3 Compile the C/C++ code

    When creating the environment at Step 1, the default machine in CompileC.sh will be set to the intended one. Also, CompileC.sh will set the directory where it is to be the environment directory and use the local HDF5 at the first instance.

    ```bash
    #  # assuming under the environment directory
    cd OpenSBLI/apps/transitional_SBLI/
    ../../../CompileC.sh
    ```
  - Step 4 Optional setup

    After successfully creating the OpenSBLI/OPS software environment after compiling and installing everything, the script will indicate that we can set up a few more environment variables by sourcing the script **OpenSBLIEnvVar** under the specified directory (using, e.g. **source ~/tmp/OpenSBLIEnv**). It will set up the $PATH so that we can run **Generate.sh** and **CompileC.sh** without using "../../../". Also, one can call **source $PYTHON** to enter into the provided Python environment for postprocessing.

    
- Example of generating, compiling, running and postprocessing an OpenSBLI app using the 1d wave problem called "wave". Follow the below outlined process:
   - cd ~/tmp/OpenSBLIEnv/OpenSBLI/apps
   - cp -r wave wave_test
   - cd wave_test
   - For generating OPS compliant cpp code from the OpenSBLI python application script **wave.py**, either do "**../../../Generate.sh wave.py**" or first do "source ~/tmp/OpenSBLIEnv" and then just do "**Generate.sh wave.py**"
   - This will create an opensbli.cpp file and other "*.h" header files
   - Now you do "**../../../CompileC.sh**", this will do the code translation using OPS to multiple intended backends. E.g., this will create a local folder named **build** and will contain by default all the executable that could be generated on a particular machine. You should see exculatable like "OpenSBLI_seq", "OpenSBLI_mpi", etc. if you do "**ll ./build**"
   - Now you can run one of these excutables in the present directory (~/tmp/OpenSBLIEnv/OpenSBLI/apps/wave_test) itself to test if it is running or not using the serial version "./build/OpenSBLI_seq" or using mpi parallel executable "mpirun -np 4 ./build/OpenSBLI_mpi".
   - This should generate an hdf5 file named "**opensbli_output.h5**". Now, to postprocess and plot the results of this file, you can run the command **source $PYTHON** to activate the Python environment. You must ensure that you have run **source ~/tmp/OpenSBLIEnv** command before doing **source $PYTHON**. Once you are inside the environment, run "**python plot.py**". This will create a folder called "simulation_plots" and will contain files "phi.pdf" and "phi_error.pdf".
   - to open these PDF files, you will need a PDF viewer, which can be installed using "**sudo apt install evince**". Then you can do "**evince simulation_plots/phi.pdf**" to visualise the line plot.
   - Once finished with postprocessing, you can now exit the python environment by doing "**conda deactivate**".


- Generate the submission script for supercomputers

  The capability is implemented by using Python 3. It requires two packages which are not shipped with Python 3 by default, i.e., Jinja2 and commentjson. The two packages are installed by default for the Python package installed with InstallPython.sh. The script will automatically install them if they are not detected.

  To use the script, we need to specify the job in a JSON format file, say "JobSubmission.json"

  ```json
  {
    "Machine":"ARCHER2",
    "Nodes" : 64,
    "TasksPerNode":128,
    "CpusPerTask":1,
    "JobName":"TGV",
    "Account":"c01-eng",
    "QoS":"standard",
    "Partition": "standard",
    "Time": "00:10:00",#HH:MM:SS
    "Program":"./TGV3DMpi Config=tgv3d1024.json"
  }
  ```
  and then run the script to generate the submission script accordingly

  ```bash
  Submit.py JobSubmission.json
  ```
- Tips

  - There is an explanation on usage in each script, which can be shown by, e.g.,

    ```bash
      ./InstallCMake.sh -h
    ```
