#!/bin/bash
##@brief Create  a working environment for OpenSBLI!
##@author Jianping Meng
##@contributors
##@details

## TODO: compiler support

function usage {
    echo "This script will download, compile and install the OPS library to a specified directory!"
    echo "./$(basename $0) -h -> Showing usage"
    echo "./$(basename $0) -b -> Specifying OpenSBLI branch (e.g.,cpc_release)"
    echo "./$(basename $0) -o -> Specifying OPS branch (e.g.,feature/HDF5Slice) "
    echo "./$(basename $0) -d -> Specifying the directory (absolute path) for creating the environment"
    echo "./$(basename $0) -p -> Use Python2"
    echo "./$(basename $0) -H -> if installing a local HDF5 library"
    echo "./$(basename $0) -c -> Specifying the compiler"
    echo "./$(basename $0) -m -> Specifying the machine type"
    echo "Machine type can be: Ubuntu ARCHER2 IRIDIS5 Fedora DAaaS"
    echo "If without specifying the machine, the script will assume all dependencies prepared!"
}
# DAaaS is a STFC platform for training, using Rocky Linux at this moment

optstring="hb:o:d:pHC:m:"
Compiler="Gnu"
Branch="version3_cpc"
Dir="$HOME/OpenSBLI"
Machine="None"
LocalHDF5="OFF"
PythonVer=""
OpsBranch="develop"

while getopts ${optstring} options; do
    case ${options} in
        h)
            usage
            exit 0
        ;;
        c)
            Compiler=${OPTARG}
        ;;
        m)
            Machine=${OPTARG}
        ;;
        d)
            Dir=${OPTARG}
        ;;
        b)
            Branch=${OPTARG}
        ;;
        o)
            OpsBranch=${OPTARG}
        ;;
        H)
            LocalHDF5="ON"
        ;;
        p)
            PythonVer="-p"
        ;;
        :)
            echo "$0: Must supply an argument to -$OPTARG." >&2
            exit 1
        ;;
        ?)
            echo "Invalid option: -${OPTARG}."
            exit 2
        ;;
    esac
done

## Install Gnu compiler provided by system
if [ $Machine == "Ubuntu" ]
then
    HDF5Lib="libhdf5-openmpi-dev"
    if [[ "${LocalHDF5}"  == "ON" ]]
    then
        HDF5Lib=""
    fi
    sudo apt install libhdf5-openmpi-dev build-essential
fi
if [ $Machine == "Fedora" ]
then
    HDF5Lib="hdf5-openmpi-devel hdf5-devel"
    if [[ "${LocalHDF5}"  == "ON" ]]
    then
        HDF5Lib=""
    fi
    sudo dnf install make automake gcc gcc-c++ kernel-devel ${HDF5Lib}
fi

if [ -d "$Dir" ]
then
   echo "The specified directory is exising, please delete it first or change to another one!"
   exit 1
fi
mkdir -p $Dir
cd $Dir ## needed by normal version
## Download the scripts

# wget -c https://github.com/jpmeng/aosh/archive/refs/heads/main.zip
wget -c https://github.com/opensbli/environment_and_run_scripts/archive/refs/heads/main.zip
unzip main.zip
rm main.zip
mv environment_and_run_scripts-main/*  .
rm -r -f environment_and_run_scripts-main
# Set the default machine in CompileC.sh changed from "Ubuntu" to "None"
sed -i "s/Machine=\"None\"/Machine=\"${Machine}\"/g" CompileC.sh
## Install local HDF5 if needed
WithHDF5=""
if [[ "${LocalHDF5}"  == "ON" ]]
then
    if [ -d "HDF5" ]
    then
        rm -r -f HDF5
    fi
    mkdir ${Dir}/HDF5
    WithHDF5="-H ${Dir}/HDF5"
    if [[ "${Machine}" == "Ubuntu" ]] ||  [[ "${Machine}" == "Fedora" ]] || [[ "${Machine}" == "DAaaS" ]]
    then
        export CC=mpicc
    fi
    ./InstallHDF5.sh -d "${Dir}/HDF5" -m ${Machine}
fi
# Python
./InstallPython.sh -d $Dir/Python ${PythonVer}
# OPS
./InstallOPS.sh -d $Dir/OPS-INSTALL -m ${Machine} ${WithHDF5} -o ${OpsBranch} -A ${Dir}
# OpenSBLI
wget -c https://github.com/opensbli/opensbli/archive/refs/heads/${Branch}.zip
FileName="$(basename -- $Branch)"
unzip "${FileName}.zip"
rm -r -f "${FileName}.zip"
DefaultOpenSBLIDir=opensbli-`echo ${Branch} | sed  's/\//-/g'`
mv ${DefaultOpenSBLIDir} ${Dir}/OpenSBLI
#wget -c https://github.com/opensbli/opensbli/archive/refs/heads/cpc_release.zip
#unzip cpc_release.zip
#rm cpc_release.zip

ScriptPath="${BASH_SOURCE:-$0}"
AbsolutScriptPath="$(realpath "${ScriptPath}")"
EnvDir="$(dirname "${AbsolutScriptPath}")"
EnvFile="OpenSBLIEnvVar"
if [[ -f ${EnvFile} ]]
then
   rm ${EnvFile}
fi
echo "To set up a few environment variables, please use source ${EnvDir}/${EnvFile}"
echo "export PYTHON=\"${EnvDir}/Python/bin/activate ${EnvDir}/Python\"" > ${EnvFile}
echo "export PATH=\$PATH:${EnvDir}" >> ${EnvFile}
echo "echo \"To use the pre-installed Python, run source \\\$PYTHON\"" >> ${EnvFile}
