!> @file
!! @brief Interpolate MPAS grid to a target grid.
!!
!! @author George Gayno NOAA/EMC

!> Initialize an FV3 model run.
!!
!! Interpolate data from the unstructured MPAS mesh to a structured grid as defined
!! by a WRF input/history file.
!!
!! This file reads a configuration namelist.
!!
!! Link the configuration namelist to ./fort.41. Then run the program
!! with preferably a number of MPI tasks which the number of MPAS cells (NOT nodes) is 
!! evenly divisible by (e.g., 36 (or 360) tasks for 36,000 input MPAS cells).
!!
!! @note For variable names “input” refers to the MPAS data input to the
!! program . “Target” refers to the  target grid.
!!
!! @author Larissa Reames CIWRO/NOAA/NSSL/FRDD
!! @return 0 for success, error code otherwise.
 program mpassit

 use mpi
 use esmf
 use utils_mod
 !use interp, only              : interp_driver

 use program_setup, only       : read_setup_namelist, LogType

 use model_grid, only          : define_target_grid,  &
                                 define_input_grid, &
                                 cleanup_input_target_grid_data

 use input_data, only         :  read_input_data
 
 use interp, only             :  interp_data
 
 use write_data, only         : write_to_file

 implicit none

 integer                      :: ierr, localpet, npets, unum, lenstr, istatus, comm
 character(100)               :: tmpstr
 logical                      :: fexist

 type(esmf_vm)                :: vm

!-------------------------------------------------------------------------
! Get namelist name from standard input
!-------------------------------------------------------------------------
 unum = COMMAND_ARGUMENT_COUNT()
 IF (unum > 0) THEN
   CALL GET_COMMAND_ARGUMENT(1, tmpstr, lenstr, istatus )
   INQUIRE(FILE=TRIM(tmpstr),EXIST=fexist)
   IF (.NOT. fexist) THEN
       call error_handler('namelist file - '//TRIM(tmpstr)//' does not exist.', -1)
   END IF
 ELSE
   print*, ' no namelist entry provided at execution, defaulting to using fort.41'
   INQUIRE(FILE='fort.41',EXIST=fexist)
   IF (.NOT. fexist) THEN
       call error_handler('namelist file fort.41 does not exist.', -1)
   END IF
 END IF

!-------------------------------------------------------------------------
! Initialize mpi
!-------------------------------------------------------------------------

 call mpi_init(ierr)

!-------------------------------------------------------------------------
! Read program configuration namelist.
!-------------------------------------------------------------------------
!
! call read_setup_namelist(filename=tmpstr)
!
!-------------------------------------------------------------------------
! Initialize mpi
!-------------------------------------------------------------------------

!  print*,"- INITIALIZE ESMF"
 call ESMF_Initialize(rc=ierr, logkindflag=LogType)
 if(ESMF_logFoundError(rcToCheck=ierr, msg=ESMF_LOGERR_PASSTHRU,  line=__LINE__,file=__FILE__)) &
    call error_handler("INITIALIZING ESMF", ierr)

 !print*,"- CALL VMGetGlobal"
 call ESMF_VMGetGlobal(vm, rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
    call error_handler("IN VMGetGlobal", ierr)

 !if (localpet==0) print*,"- CALL VMGet"
 call ESMF_VMGet(vm, localPet=localpet, petCount=npets, mpiCommunicator=comm, rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
    call error_handler("IN VMGet", ierr)

 if (localpet==0) print*,'- NPETS IS  ',npets
 !print*,'- LOCAL PET ',localpet

!-------------------------------------------------------------------------
! Read program configuration namelist.
!-------------------------------------------------------------------------

 call read_setup_namelist(filename=tmpstr)

!-------------------------------------------------------------------------
! Create esmf grid objects for input and target grids.
!-------------------------------------------------------------------------

 call define_target_grid(localpet, npets)
 
 call define_input_grid(localpet, npets)

!-------------------------------------------------------------------------
! Read data from input file
!-------------------------------------------------------------------------

call read_input_data(localpet, npets, comm)

!-------------------------------------------------------------------------
! Interpolate fields
!-------------------------------------------------------------------------

 call interp_data(localpet)
 
!-------------------------------------------------------------------------
! Write data to file
!-------------------------------------------------------------------------

 call write_to_file(localpet)
 
!-------------------------------------------------------------------------
! Finish up
!-------------------------------------------------------------------------

 call cleanup_input_target_grid_data(localpet)

 if (localpet==0) print*,"- CALL ESMF_finalize"
 call ESMF_finalize(endflag=ESMF_END_KEEPMPI, rc=ierr)

 call mpi_finalize(ierr)

 print*,"- DONE."

 end program mpassit
