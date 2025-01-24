help([[
This module loads libraries for MPASSIT
]])

whatis([===[Loads libraries for rrfs-workflow ]===])
prepend_path("MODULEPATH", "/ncrc/proj/epic/spack-stack/c6/spack-stack-1.6.0/envs/unified-env/install/modulefiles/Core")
load("stack-intel/2023.2.0")
load("stack-cray-mpich/8.1.29")
load("parallel-netcdf/1.12.2")
load("parallelio/2.5.10")

load("jasper/2.0.32")
load("zlib/1.2.13")
load("libpng/1.6.37")
load("esmf/8.5.0")

setenv("CMAKE_C_COMPILER", "cc")
setenv("CMAKE_CXX_COMPILER", "CC")
setenv("CMAKE_Fortran_COMPILER", "fn")
