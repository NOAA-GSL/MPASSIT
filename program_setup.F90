!> @file
!! @brief Set up program execution.
!! @author Larissa Reames CIWRO/NOAA/NSSL/FRDD

!> This module contains code to read the setup namelist file.
!!
!! @author Larissa Reames CIWRO/NOAA/NSSL/FRDD
 module program_setup

 use esmf
 use ESMF_LogPublicMod
 use constants_module
 use utils_mod
 use misc_definitions_module

 implicit none

 private

 type(ESMF_LogKind_Flag), public :: LogType 
 
 ! Namelist variables
 character(len=500), public      :: grid_file_input_grid = "NULL" !< Full path of MPAS file containing grid information
 character(len=500), public      :: diag_file_input_grid = "NULL" !< Full path of input diagnostic MPAS data
 character(len=500), public      :: hist_file_input_grid = "NULL" !< Full path of input history MPAS data
 character(len=500), public      :: file_target_grid = "NULL"     !<Full path of file containing target 
                                                                  !<grid information for target_grid_type='file'
 character(len=500), public      :: output_file = "NULL"          !< Full path of output file
 
 logical, public                 :: interp_diag = .false. !< Read data from diag file?
 logical, public                 :: interp_hist = .false. !< Read data from hist file?
 logical, public                 :: wrf_mod_vars = .false.!< Whether to modify variable values/dimensions 
                                                          !< to conform to WRF format. Set to true for
                                                          !< UPP-compatible output
 character(len=500), public      :: target_grid_type      !< Grid type to interpolate data to
                                                          !< Valid options: 'file', 'lambert',
                                                          !< 'mercator','polar',lat-lon'     
 character(len=500), public      :: block_decomp_file = "NULL"  !< Full path to MPAS grid-specific block decomposition file
 real(kind=4), public      :: missing_value = -9999.   !< Value to be used for mapped data outside of input data

 !! These entries are only valid when target_grid_type is not 'file'
 logical, public                 :: is_regional = .true.  !< Is the output grid regional or global? 
                                                          !< Default: True   
 integer, public                 :: i_target              !< # staggered east-west grid points in target grid
 integer, public                 :: j_target              !< # staggered north-south grid points in target grid                                   
 real, public                    :: truelat1 = NAN        !< First true latitude (all projections)
 real, public                    :: truelat2 = NAN        !< Second true latitude (LCC only)
 real, public                    :: stand_lon = NAN       !< Longitude parallel to y-axis (-180->180E)
 real, public                    :: dx = NAN              !< Grid cell east-west dimension(meters or deg for 
                                                          !< target_grid_type='lat-lon')
 real, public                    :: dy = NAN              !< Grid cell north-south dimension(meters or deg for 
                                                          !< target_grid_type='lat-lon')
 real, public                    :: ref_lat               !< Latitude of reference point
 real, public                    :: ref_lon               !< Longitude of reference point
 real, public                    :: ref_x                 !< Grid-relative e-w index of reference point
                                                          !< Defaults to grid center (nx/2)
 real, public                    :: ref_y                 !< Grid-relative n-s index of reference point
                                                          !< Defaults to grid center (ny/2) 
 real, public                    :: pole_lat              !< Latitude of pole for target grid projection
 real, public                    :: pole_lon              !< Longitude of pole for target grid projection 
 
 !These aren't namelist variables but they're created directly from them
 real, public                    :: dxkm                  !< grid-cell east-west dimension (meters)
 real, public                    :: dykm                  !< grid-cell north-south dimension (meters)
 real, public                    :: dlondeg               !< grid-cell east-west dimension (deg)
 real, public                    :: dlatdeg               !< grid-cell north-south dimension (deg)
 real, public                    :: known_lat             !< Latitude of reference point
 real, public                    :: known_lon             !< Longitude of reference point
 real, public                    :: known_x               !< Grid-relative e-w index of reference point
 real, public                    :: known_y               !< Grid-relative n-s index of reference point 
 integer, public                 :: proj_code             !< Integer code corresponding to the requested
                                                          !< target grid projection type
 character(len=500), public      :: map_proj_char         !< Map projection name
 logical, public                 :: interp_as_bundle = .true. !< If true, use ESMF_FieldBundleRegridStore and interpolate the fields as a bundle
                                                              !< If false, use ESMFRegridStore and interpolate individual fields.
                                                              !< .false. seems faster and less memory intensive
                                                              !< Currently, only applies to conservative regridding
                                                           

 public :: read_setup_namelist

 contains

!> Reads program configuration namelist.
!!
!! @param filename the name of the configuration file (defaults to
!! ./fort.41).
!! @author Larissa Reames CIWRO/NOAA/NSSL/FRDD
 subroutine read_setup_namelist(unum, filename)
 implicit none

 character(len=*), intent(in), optional :: filename
 integer, intent(in), optional          :: unum
 character(:), allocatable              :: filename_to_use
 character(len=500)                     :: map_proj
 integer                                :: unit_to_use
 logical                                :: esmf_log,decomp_exists

 integer                                :: is, ie, ierr
 
 !Namelist variables that are used to create global variables
 integer                                :: nx,ny

 namelist /config/ grid_file_input_grid, diag_file_input_grid, hist_file_input_grid, &
            file_target_grid, output_file, interp_diag, interp_hist, &
            wrf_mod_vars, esmf_log,target_grid_type,nx,ny,dx,dy,ref_lat,ref_lon,ref_x,ref_y,&
            truelat1,truelat2,stand_lon,is_regional,pole_lat,pole_lon, interp_as_bundle,block_decomp_file, &
            missing_value

  ref_x = NAN
  ref_y = NAN
  ref_lat = NAN
  ref_lon = NAN
  dx = NAN
  dy = NAN
  pole_lat = 90.0
  pole_lon = 0.0
  nx = 0
  ny = 0                   

 !print*,"- READ SETUP NAMELIST"

 if (present(filename)) then
    filename_to_use = filename
 else
    filename_to_use = "./fort.41"
 endif

 if (present(unum)) then
     unit_to_use = unum
 else
     unit_to_use = 41
 endif

 open(unit_to_use, file=filename_to_use, iostat=ierr)
 if (ierr /= 0) call error_handler("OPENING SETUP NAMELIST.", ierr)
 read(unit_to_use, nml=config, iostat=ierr)
 if (ierr /= 0) call error_handler("READING SETUP NAMELIST.", ierr)
 close (unit_to_use)
 
 if (esmf_log) then
   LogType = ESMF_LOGKIND_MULTI_ON_ERROR
 else
   LogType = ESMF_LOGKIND_NONE
 endif 
 
 !if (block_decomp_file=='NULL') then
   !call error_handler("block_decomp_file IS REQUIRED BUT IS MISSING IN NAMELIST.", -1)
 !else
 if(.not. block_decomp_file=='NULL') then
   inquire(file=block_decomp_file, exist=decomp_exists)
   if (.not. decomp_exists) then
      call error_handler("block_decomp_file DOES NOT EXIST.",-1)
   endif
 endif
 
 if (trim(target_grid_type) .ne. 'file') then
   dxkm = dx
   dykm = dy
 
   known_lat = ref_lat
   known_lon = ref_lon
   known_x = ref_x
   known_y = ref_y
   i_target = nx-1
   j_target = ny-1
 
   map_proj = to_upper(target_grid_type)
   !print*, map_proj 
   ! Assign parameters to module variables
   if ((index(map_proj, 'LAMBERT') /= 0) .and. &
      (len_trim(map_proj) == len('LAMBERT'))) then
     proj_code = PROJ_LC 
     map_proj_char = 'Lambert Conformal'

   else if ((index(map_proj, 'MERCATOR') /= 0) .and. &
           (len_trim(map_proj) == len('MERCATOR'))) then
     proj_code = PROJ_MERC 
     map_proj_char = 'Mercator'

   else if ((index(map_proj, 'POLAR') /= 0) .and. &
           (len_trim(map_proj) == len('POLAR'))) then
     proj_code = PROJ_PS 
     map_proj_char = "Polar Stereographic"

   else if ((index(map_proj, 'LAT-LON') /= 0) .and. &
           (len_trim(map_proj) == len('LAT-LON'))) then
     proj_code = PROJ_LATLON
     map_proj_char = 'Lat/Lon'
   else
     call error_handler('In namelist, invalid target_grid_type specified. Valid '// &
                  'projections are "lambert", "mercator", "polar", and '// &
                  '"lat-lon".',ERROR_CODE)
   end if
 
 
   if (proj_code == PROJ_LATLON) then
     ! If no dx,dy specified, assume global grid
     if (dx == NAN .and. dy == NAN) then
        if (is_regional) then
            call error_handler("For lat-lon projection, if dx/dy are not specified "// &
            "a global grid is assumed. Please set dx/dy if a regional grid is desired, "//&
            "or change is_regional to .false. if a global grid is desired.", ERROR_CODE)
        endif
        dlondeg = 360. / (i_target)   ! Here, we really do not want e_we-s_we+1
        dlatdeg = 180. / (j_target)   ! Here, we really do not want e_we-s_we+1
        known_x = 1.
        known_y = 1.
        known_lon = stand_lon + dlondeg/2.
        known_lat = -90. + dlatdeg/2.
        dxkm = EARTH_RADIUS_M * PI * 2.0 / (i_target)
        dykm = EARTH_RADIUS_M * PI       / (j_target)

     ! If dx,dy specified, however, assume regional grid
     else
        if (.not. is_regional) then
            call error_handler("For lat-lon projection, if dx/dy are specified "// &
            "a regional grid is assumed. Please unset dx/dy if a global grid is desired, "//&
            "or change is_regional to .true. if a regional grid is desired.", ERROR_CODE)
        endif
        dlatdeg = dy
        dlondeg = dx
        dxkm = dlondeg * EARTH_RADIUS_M * PI * 2.0 / 360.0
        dykm = dlatdeg * EARTH_RADIUS_M * PI * 2.0 / 360.0
        !print*, "dxkm, dykm = ", dxkm, dykm
        if (known_lat == NAN .or. known_lon == NAN) then
           call error_handler('For lat-lon projection, if dx/dy are specified, '// &
                    'a regional domain is assumed, and a ref_lat,ref_lon must also be specified',ERROR_CODE)
        end if
     end if
   end if
 
 ! Manually set truelat2 = truelat1 if truelat2 not specified for Lambert
  if (proj_code==PROJ_LC .and. truelat2 == NAN) then
     if (truelat1 == NAN) call error_handler("No TRUELAT1 specified for Lambert conformal projection.",ERROR_CODE) 
     truelat2 = truelat1
  endif
  
  ! If the user hasn't supplied a known_x and known_y, assume the center of domain 1
  if (known_x == NAN .and. known_y == NAN) then
    known_x = real(i_target+1) / 2.
    known_y = real(j_target+1) / 2.
   ! print*, known_x, known_y
  else if (known_x == NAN .or. known_y == NAN) then
    call error_handler('In namelist, neither or both of ref_x, ref_y must be specified.',ERROR_CODE)
  end if 
 endif

 return
 
 end subroutine read_setup_namelist

 end module program_setup
