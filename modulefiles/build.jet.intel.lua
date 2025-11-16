help([[
This module loads libraries for MPASSIT
]])

whatis([===[Loads libraries for MPASSIT ]===])
prepend_path("MODULEPATH","/mnt/lfs6/BMC/nrtrr/spack-stack/1.9.2/envs/ue-oneapi-2024.2.1/install/modulefiles/Core")

load("stack-oneapi/2024.2.1")
load("stack-intel-oneapi-mpi/2021.13")

load("cmake/3.27.9")
load("esmf/8.8.0")

setenv("CMAKE_C_COMPILER", "mpiicc")
setenv("CMAKE_CXX_COMPILER", "mpiicpc")
setenv("CMAKE_Fortran_COMPILER", "mpiifort")
