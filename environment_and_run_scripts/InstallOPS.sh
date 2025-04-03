#!/bin/bash
##@brief Download, compile and install OPS library
##@author Jianping Meng
##@contributors Pushpender Sharma Teja Ala
##@details

function usage {
    echo "This script will download, compile and install the OPS library to a specified directory!"
    echo "./$(basename $0) -h -> Showing usage"
    echo "./$(basename $0) -d -> Specifying the directory for installation"
    echo "./$(basename $0) -c -> Specifying the compiler"
    echo "./$(basename $0) -H -> Specifying the HDF5 directory"
    echo "./$(basename $0) -o -> Specifying the branch (e.g.,feature/HDF5Slice)"
    echo "./$(basename $0) -m -> Specifying the machine type"
    echo "./$(basename $0) -A -> Copy the CMake file for compiling applications to the directory"
    echo "Machine type can be: Ubuntu ARCHER2 IRIDIS5 Fedora DAaaS"
    echo "If without specifying the machine, the script will assume all dependencies prepared!"
}
optstring="hc:m:A:d:H:o:"
Compiler="Gnu"
Dir="$HOME/OPS_INSTALL"
Machine="None"
HDF5Root=""
Branch="develop"
AppCMakeDir=""

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
        A)
            AppCMakeDir=${OPTARG}
        ;;
        d)
            Dir=${OPTARG}
        ;;
        H)
            HDF5Root="-DHDF5_ROOT=${OPTARG}"
        ;;
        o)
            Branch="${OPTARG}"
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

if [ $# -eq 0 ]
then
    echo "This script will download, compile with ${Compiler} and install the OPS library to to ${Dir}!"
fi

# OPTIMISATION=-DCFLAG="-ftree-vectorize -funroll-loops"
# echo ${OPTIMISATION}
# #-DCXXFLAG="-ftree-vectorize -funroll-loops"

OPTIMISATION="-DCFLAG=-ftree-vectorize -DCXXFLAG=-ftree-vectorize"

if ! grep "CreateOpenSBLIEnv" /proc/$PPID/cmdline
then
    if [ $Machine == "Ubuntu" ]
    then
        HDF5Lib="libhdf5-openmpi-dev"
        if [ -z ${HDF5Root} ]
        then
            HDF5Lib=""
        fi
        sudo apt install libhdf5-mpi-dev build-essential ${HDF5Lib}
    fi
    if [ $Machine == "Fedora" ]
    then
        HDF5Lib="hdf5-openmpi-devel hdf5-devel"
        if [ -z ${HDF5Root} ]
        then
            HDF5Lib=""
        fi
        sudo dnf install make automake gcc gcc-c++ kernel-devel ${HDF5Lib}
    fi
fi

if [ $Machine == "CIRRUS" ]
then
    module load nvidia/nvhpc/22.2
    module load cmake/3.22.1
    OPTIMISATION=""
    if [ -z ${HDF5Root} ]
    then
        module load hdf5parallel/1.12.0-nvhpc-openmpi
    fi
fi

if [ $Machine == "ARCHER2" ]
then
    module purge PrgEnv-cray
    module load load-epcc-module
    module load cmake/3.21.3
    module load PrgEnv-gnu
    module load cray-python/3.9.13.1

    if [ -z ${HDF5Root} ]
    then
        module load cray-hdf5-parallel
    fi
fi

if [ $Machine == "IRIDIS5" ]
then
    module load gcc/6.4.0
    if [ -z ${HDF5Root} ]
    then
        module load hdf5/1.10.2/gcc/parallel
    fi
    module load cuda/10.0
    module load cmake
    module load python/3.9.7
fi

#https://github.com/OP-DSL/OPS/archive/refs/heads/feature/HDF5Slice.zip
#https://github.com/OP-DSL/OPS/archive/refs/heads/BurgerEquation.zip
#OPS-feature-HDF5Slice/
#OPS-BurgerEquation

wget -c https://github.com/OP-DSL/OPS/archive/refs/heads/${Branch}.zip
FileName="$(basename -- $Branch)"
unzip "${FileName}.zip"
rm -r -f "${FileName}.zip"
SourceDir=OPS-`echo ${Branch} | sed  's/\//-/g'`
cd ${SourceDir}
if [ ! -z ${AppCMakeDir} ]
then
    cp apps/c/CMakeLists.txt ${AppCMakeDir}
fi
mkdir build
cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=${Dir} -DCMAKE_BUILD_TYPE=Release ${OPTIMISATION} -DBUILD_OPS_APPS=OFF ${HDF5Root}
cmake --build . -j 4
cmake --install .
cd ../../
rm -r -f ${SourceDir}
