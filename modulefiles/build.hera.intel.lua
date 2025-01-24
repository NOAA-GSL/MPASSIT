help([[
This module loads libraries for MPASSIT
]])

whatis([===[Loads libraries for rrfs-workflow ]===])
prepend_path("MODULEPATH", "/scratch1/NCEPDEV/nems/role.epic/spack-stack/spack-stack-1.6.0/envs/unified-env-rocky8/install/modulefiles/Core")

load("stack-intel/2021.5.0")
load("stack-intel-oneapi-mpi/2021.5.1")
load("cmake/3.23.1")
load("gnu")
load("intel/2022.1.2")
load("impi/2022.1.2")

load("pnetcdf/1.7.0")
load("szip")
load("hdf5parallel/1.10.6")
load("netcdf-hdf5parallel/4.7.4")
load("jasper/2.0.32")
load("zlib/1.2.13")
load("libpng/1.6.37")
load("esmf/8.5.0")

setenv("CMAKE_C_COMPILER", "mpiicc")
setenv("CMAKE_CXX_COMPILER", "mpiicpc")
setenv("CMAKE_Fortran_COMPILER", "mpiifort")
