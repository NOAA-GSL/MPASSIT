help([[
This module loads libraries for MPASSIT
]])

whatis([===[Loads libraries for MPASSIT ]===])
prepend_path("MODULEPATH", "/apps/contrib/spack-stack/spack-stack-1.9.3/envs/ue-oneapi-2024.2.1/install/modulefiles/Core")

load("stack-oneapi/2024.2.1")
load("stack-intel-oneapi-mpi/2021.13")

load("cmake/3.27.9")
load("esmf/8.8.0")

setenv("CMAKE_C_COMPILER", "mpiicc")
setenv("CMAKE_CXX_COMPILER", "mpiicpc")
setenv("CMAKE_Fortran_COMPILER", "mpiifort")
