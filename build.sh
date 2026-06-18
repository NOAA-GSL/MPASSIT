#! /usr/bin/env bash
#
# Author: Larissa Reames CIWRO/NOAA/NSSL/FRDD

#set -eux

target=${1:-""}
compiler=${2:-"intel"}
debug=${3:-"true"}

# If target is not set
if [[ "$target" == "" ]]; then
    source ./detect_machine.sh
    target=${MACHINE}
fi

echo "target=$target, compiler=$compiler"

if [[ "$target" == "" ]]; then
  echo "target is not set and the platform name cannot be detected automatically"
  exit 1
fi

# Check for platform/compiler configuration file
if [[ ! -f modulefiles/build.$target && ! -f modulefiles/build.$target.$compiler.lua && ! -f modulefiles/build.$target.$compiler ]]; then
    echo "Platform ${target} configuration file not found in ./modulefiles, neither build.$target nor build.$target.$compiler.lua"
    exit 1
fi

if [[ "$target" == "vecna" ]]; then
    echo "Use platform configuration file: build.$target.$compiler"
    source ./modulefiles/build.$target.$compiler > /dev/null
elif [[ "$target" == "linux.*" || "$target" == "macosx.*" ]]; then
    unset -f module
    echo "Use platform configuration file: build.$target"
    source ./modulefiles/build.$target > /dev/null
else
    echo "Use platform configuration file: build.$target.$compiler.lua"
    module use ./modulefiles
    module load build.$target.$compiler.lua
    module list
fi

CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=../ -DEMC_EXEC_DIR=ON -DBUILD_TESTING=OFF"
if [[ "$target" == "wcoss2" ]]; then
    CMAKE_FLAGS="${CMAKE_FLAGS} -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC -DCMAKE_Fortran_COMPILER=ftn"
elif [[ "$compiler" == "intel-llvm" ]]; then
    CMAKE_FLAGS="${CMAKE_FLAGS} -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx -DCMAKE_Fortran_COMPILER=ifx"
fi

if [[ "${debug}" == "true" ]]; then
    CMAKE_FLAGS="${CMAKE_FLAGS} -DCMAKE_BUILD_TYPE=Debug"
else
    CMAKE_FLAGS="${CMAKE_FLAGS} -DCMAKE_BUILD_TYPE=Release"
fi

# for a clean build folder
rm -fr ./build
mkdir ./build && cd ./build || exit 0

# do the building
cmake .. ${CMAKE_FLAGS}

make -j 8 VERBOSE=1
make install

exit 0
