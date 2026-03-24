!> @file
!! @brief Read atmospheric and surface data from MPAS diag, hist, or forecast files.
!! @author  Larissa Reames CIWRO/NOAA/NSSL

!> Read atmospheric and/or surface input MPAS grid.
!! Supported formats include INITIalization file, forecast diag file (only 2d diagnostic
!! fields), and full 3d forecast files.
!!
!! Public variables are defined below: "input" indicates field
!! associated with the input grid.
!!
!! @author  Larissa Reames CIWRO/NOAA/NSSL

 module input_data

 use esmf
 use netcdf
 use utils_mod
 use program_setup, only          : hist_file_input_grid, &
                                    diag_file_input_grid, &
                                    grid_file_input_grid, &
                                    interp_diag, interp_hist, &
                                    wrf_mod_vars

 use model_grid, only             : input_grid,        &
                                    nCells_input, nVert_input,  &
                                    nz_input, nzp1_input, &
                                    strlen, &
                                    valid_time, config_dt, &
                                    start_time, lsm_scheme, &
                                    mp_scheme, conv_scheme, &
                                    cell_latitude_input_grid, &
                                    cell_longitude_input_grid, &
                                    zgrid_input_grid, &
                                    input_diag_bundle, &
                                    target_diag_names, &
                                    input_hist_bundle_2d_cons, &
                                    input_hist_bundle_2d_patch, &
                                    input_hist_bundle_2d_nstd, &
                                    input_hist_bundle_3d_nz, &
                                    input_hist_bundle_3d_nzp1, &
                                    input_hist_bundle_3d_vert, &
                                    input_hist_bundle_soil, &
                                    target_hist_names_2d_cons, &
                                    target_hist_names_2d_nstd, &
                                    target_hist_names_2d_patch, &
                                    target_hist_names_3d_nz, &
                                    target_hist_names_3d_nzp1, &
                                    target_hist_names_3d_vert, &
                                    target_hist_names_soil, &
                                    target_diag_units, &
                                    target_hist_units_2d_cons, &
                                    target_hist_units_2d_nstd, &
                                    target_hist_units_2d_patch, &
                                    target_hist_units_3d_nzp1, &
                                    target_hist_units_3d_vert, &
                                    target_hist_units_3d_nz, &
                                    target_hist_units_soil, &
                                    target_diag_longname, &
                                    target_hist_longname_2d_cons, &
                                    target_hist_longname_2d_nstd, &
                                    target_hist_longname_2d_patch, &
                                    target_hist_longname_3d_nzp1, &
                                    target_hist_longname_3d_nz, &
                                    target_hist_longname_3d_vert, &
                                    target_hist_longname_soil, &
                                    n_diag_fields, &
                                    n_hist_fields_2d_cons, &
                                    n_hist_fields_2d_nstd, &
                                    n_hist_fields_2d_patch, &
                                    n_hist_fields_3d_vert, &
                                    n_hist_fields_3d_nz, &
                                    n_hist_fields_3d_nzp1, &
                                    n_hist_fields_3d_vert, &
                                    n_hist_fields_soil, &
                                    elemIDs, nCellsPerPET, &
                                    nodeIDs,diag_out_interval, &
                                    nNodesPerPET, &
                                    u_input_grid, &
                                    v_input_grid, &
                                    do_u_interp, &
                                    do_v_interp, &
                                    do_u10_interp, &
                                    do_v10_interp, &
                                    u10_ind, v10_ind
 implicit none

 private

 integer, public                        :: nsoil_input
                                           !< number of input soil levels
 real(esmf_kind_r8), allocatable , public :: zs_target_grid(:,:)
                                          !< soil center depth, target grid
 real(esmf_kind_r8), allocatable , public :: dzs_target_grid(:,:)
                                          !< soil layer thickness, target grid

 public :: read_input_data

 contains

!> Calculates mapping from NetCDF indices to ESMF indices.
!!
!! This calculates the mapping used in scatter_ncvar_3d and scatter_ncvar_2d.
!! It maps from NetCDF file indices to ESMF indices to prepare for a FieldScatter
!!
!! @param[in] localpet  ESMF local persistent execution thread
!! @param[in] npets     number of PETs
!! @param[in] comm      MPI communicator for the PETs
!! @param[out] mapping  The output mapping.
!! @author Sam Trahan NOAA/GSL/CIRES
 subroutine gather_global_cell_mapping(localpet, npets, comm, mapping)
   use mpi
   implicit none
   integer, intent(in) :: localpet, npets
   integer, intent(inout) :: mapping(nCells_input)
   integer :: comm

   integer :: mycount, counts(npets), ierr, displs(npets), i

   mycount = size(elemIDs)
   call MPI_Allgather(mycount, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, comm, ierr)

   displs(1) = 0
   do i=2,npets
      displs(i) = displs(i-1) + counts(i-1)
   enddo

   call MPI_Allgatherv(elemIDs, mycount, MPI_INTEGER, mapping, counts, displs, MPI_INTEGER, comm, ierr)
 end subroutine gather_global_cell_mapping


!> Reads a 3D variable on rank 0 and scatters (ESMF_FieldScatter) to all ranks.
!!
!! @param[in] field  The ESMF field to read and scatter.
!! @param[in] localpet  ESMF local persistent execution thread
!! @param[inout] read_buffer  A buffer for reading from the NetCDF. Must match the shape in the file.
!! @param[inout] scatter_buffer  Input to ESMF_FieldScatter. Must match the shape of the ESMF variable for all PETs.
!! @param[in] mapping  The mapping from gather_global_cell_mapping.
!! @param[in] ncid  From NetCDF, the ncid of the file to read. Only significant on rank 0.
!! @param[in] id_var  From NetCDF, the variable id of the variable to read. Only significant on rank 0.
!! @param[in] nk  Extent of the vertical dimension.
!! @author Sam Trahan NOAA/GSL/CIRES
 subroutine scatter_ncvar_3d(field, localpet, read_buffer, scatter_buffer, mapping, ncid, id_var, nk)
   implicit none
   type(esmf_field) :: field
   integer, intent(in) :: localpet, ncid, id_var, mapping(:), nk
   real(esmf_kind_r8) :: read_buffer(:,:,:)
   real(esmf_kind_r8) :: scatter_buffer(:,:)

   integer :: j, k, error, rc

   if (localpet==0) then
      error=nf90_get_var(ncid, id_var, read_buffer)
      call netcdf_err(error, 'reading field' )
      do k=1,nk
         do j=1,nCells_input
            scatter_buffer(j,k) = read_buffer(k,mapping(j),1)
         enddo
      enddo
   endif
   call ESMF_FieldScatter(field, scatter_buffer, rootPet=0, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldScatter 3D", rc)
 end subroutine scatter_ncvar_3d

!> Reads a 2D variable on rank 0 and scatters (ESMF_FieldScatter) to all ranks.
!!
!! @param[in] field  The ESMF field to read and scatter.
!! @param[in] localpet  ESMF local persistent execution thread
!! @param[inout] read_buffer  A buffer for reading from the NetCDF. Must match the shape in the file.
!! @param[inout] scatter_buffer  Input to ESMF_FieldScatter. Must match the shape of the ESMF variable for all PETs.
!! @param[in] mapping  The mapping from gather_global_cell_mapping.
!! @param[in] ncid  From NetCDF, the ncid of the file to read. Only significant on rank 0.
!! @param[in] id_var  From NetCDF, the variable id of the variable to read. Only significant on rank 0.
!! @author Sam Trahan NOAA/GSL/CIRES
 subroutine scatter_ncvar_2d(field, localpet, read_buffer, scatter_buffer, mapping, ncid, id_var)
   implicit none
   type(esmf_field) :: field
   integer, intent(in) :: localpet, ncid, id_var, mapping(:)
   real(esmf_kind_r8) :: read_buffer(:)
   real(esmf_kind_r8) :: scatter_buffer(:)

   integer :: j, error, rc

   if (localpet==0) then
      error=nf90_get_var(ncid, id_var, read_buffer)
      call netcdf_err(error, 'reading field' )
      do j=1,nCells_input
         scatter_buffer(j) = read_buffer(mapping(j))
      enddo
   endif
   call ESMF_FieldScatter(field, scatter_buffer, rootPet=0, rc=rc)
   if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldScatter 2D", rc)
 end subroutine scatter_ncvar_2d

!> Read input grid data driver.
!!
!! @param[in] localpet  ESMF local persistent execution thread
!! @author Larissa Reames CIWRO/NOAA/NSSL
 subroutine read_input_data(localpet, npets, comm)
 use mpi
 implicit none

 integer, intent(in)             :: localpet, npets, comm

 integer, allocatable, target :: mapping(:)

 allocate(mapping(nCells_input))
 call gather_global_cell_mapping(localpet, npets, comm, mapping)

 if (interp_diag) then
    call read_input_diag_data(localpet, npets, mapping)
 endif

 if (interp_hist) then
    call read_input_hist_data(localpet, npets, mapping)
 endif

 if (.not. interp_diag .and. .not. interp_hist) then
    call error_handler(" SET INTERP_DIAG AND/OR INTERP_HIST TO TRUE TO OBTAIN OUTPUT", -1)
 endif

 deallocate(mapping)

 end subroutine read_input_data

!> Read input grid diag data.
!!
!! @param[in] localpet  ESMF local persistent execution thread
!! @author Larissa Reames CIWRO/NOAA/NSSL
 subroutine read_input_diag_data(localpet, npets, mapping)

 character(len=500)              :: the_file
 character(len=50)               :: vname

 integer, intent(in)             :: localpet, npets, mapping(:)
 integer                         :: error, ncid, rc
 integer                         :: id_dim
 integer                         :: id_var, i, j, k, nodes, ndims

 type(esmf_field),allocatable, target    :: fields(:)

 real(esmf_kind_r8), allocatable, target :: dummy(:), scatter(:), dummy2(:,:,:), scatter2(:,:)

 real(esmf_kind_r8), pointer     :: varptr(:), varptr2(:,:)

 call init_input_diag_fields(localpet)

 if (localpet==0) print*,"- READ INPUT DIAG DATA."


 the_file = trim(diag_file_input_grid)
 error=nf90_open(trim(the_file),nf90_nowrite,ncid)
 call netcdf_err(error, 'opening: '//trim(the_file) )


!---------------------------------------------------------------------------
! Initialize esmf atmospheric fields.
!---------------------------------------------------------------------------
 allocate(fields(n_diag_fields))
 allocate(target_diag_units(n_diag_fields))
 allocate(target_diag_longname(n_diag_fields))
 call ESMF_FieldBundleGet(input_diag_bundle, fieldList=fields, &
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldBundleGet", rc)

 call ESMF_MeshGet(input_grid, nodeCount = nodes, rc=rc)
 if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN MeshGet", rc)

 if(localpet == 0) then
    allocate(dummy(nCells_input))
    allocate(scatter(nCells_input))
    allocate(dummy2(nz_input, nCells_input,1))
    allocate(scatter2(nCells_input, nz_input))
 else
    allocate(dummy(1))
    allocate(scatter(1))
    allocate(dummy2(1,1,1))
    allocate(scatter2(1,1))
 endif

 do i = 1,n_diag_fields
    
    call ESMF_FieldGet(fields(i), name=vname, rc=rc)
    if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldGet", rc)
    if(vname=='u10') then
        do_u10_interp = 1
        u10_ind = i
    endif
    if(vname=='v10') then 
        do_v10_interp = 1
        v10_ind = i
    endif

    if (localpet==0) print*,"- READ ", trim(vname)
    error=nf90_inq_varid(ncid, trim(vname), id_var)
    call netcdf_err(error, 'reading field id - '//trim(vname) )
    error=nf90_inquire_variable(ncid, id_var, ndims=ndims)
    call netcdf_err(error, 'reading variable dims' )
    if (localpet==0) print*,"- READ INPUT FOR ", trim(vname)
    if (ndims == 2) then
      call scatter_ncvar_2d(fields(i), localpet, dummy, scatter, mapping, ncid, id_var)
    elseif (ndims == 3) then
      call scatter_ncvar_3d(fields(i), localpet, dummy2, scatter2, mapping, ncid, id_var, nz_input)
    endif
    error=nf90_get_att(ncid,id_var,'units',target_diag_units(i))
    call netcdf_err(error, 'reading field units' )
    error=nf90_get_att(ncid,id_var,'long_name',target_diag_longname(i))
    call netcdf_err(error, 'reading field long_name' )
 enddo

 if (localpet==0) print*,'- TRY TO READ GLOBAL ATTRIBUTE START TIME FROM DIAG FILE'
 error = nf90_get_att(ncid,NF90_GLOBAL,'config_start_time',start_time)
 if ( error .ne. 0 ) then
    if (interp_hist) then
       if (localpet==0) then
           print*,'GETTING GLOBAL ATTRIBUTE START TIME FROM DIAG FILE FAILED'
           print*,'...THAT IS OK...TRY TO GET FROM HISTORY FILE'
       endif
    else
       call netcdf_err(error, 'reading config_start_time')
    endif
 endif

 if (localpet==0) print*,'- TRY TO READ GLOBAL ATTRIBUTE CONFIG_DT FROM DIAG FILE'
 error = nf90_get_att(ncid,NF90_GLOBAL,'config_dt',config_dt)
 if (error .ne. 0) then
    if (localpet==0) then
        print*,'GETTING GLOBAL ATTRIBUTE CONFIG_DT FROM DIAG FILE FAILED'
        print*,'...THAT IS OK...SET TO 0 AND TRY TO GET FROM HISTORY FILE'
    endif
    config_dt = 0.0
 endif

 if (localpet==0) print*,'- read output_interval from diag file'
 error = nf90_get_att(ncid,NF90_GLOBAL,'output_interval',diag_out_interval)
 if (error .ne. 0) then
   if (localpet==0) print*, 'error reading output_interval from diag file, setting to 0'
   diag_out_interval = 0
 endif

 if (localpet==0) print*, "- GETTING xtime"
 error = nf90_inq_dimid(ncid,'StrLen',id_var)
 call netcdf_err(error, 'reading strlen dim id')
 error = nf90_inquire_dimension(ncid,id_var,len=strlen)
 allocate(valid_time(1,strlen))
 if (localpet==0) print*, "... xtime strlen = ", strlen
 error = nf90_inq_varid(ncid,'xtime',id_var)
 call netcdf_err(error, 'reading xtime id')
 error = nf90_get_var(ncid, id_var, valid_time)
 call netcdf_err(error, 'getting xtime')

 error = nf90_close(ncid)


 deallocate( dummy)

 end subroutine read_input_diag_data

 subroutine init_input_diag_fields(localpet)

    implicit none

    integer, intent(in)           :: localpet
    integer                       :: i, rc
    type(esmf_field),allocatable  :: diag_fields(:)
    character(50), allocatable    :: input_diag_names(:)
    character(50)                 :: fname

    fname = 'diaglist'

    call read_varlist(localpet,fname,n_diag_fields,input_diag_names, target_diag_names)

    allocate(diag_fields(n_diag_fields))
    if (localpet==0) print*,"- INITIALIZE INPUT DIAG FIELDS."
    do i = 1, n_diag_fields
        if (input_diag_names(i) == "refl10cm") then
          if (localpet==0) print*, "- INIT FIELD ", input_diag_names(i)
          diag_fields(i) = ESMF_FieldCreate(input_grid, &
                            typekind=ESMF_TYPEKIND_R8, &
                            meshloc=ESMF_MESHLOC_ELEMENT, &
                            ungriddedLBound=(/1/), &
                            ungriddedUBound=(/nz_input/), &
                            name=input_diag_names(i), rc=rc)
          if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__))&
          call error_handler("IN FieldCreate", rc)
        else
          if (localpet==0) print*, "- INIT FIELD ", input_diag_names(i)
          diag_fields(i) = ESMF_FieldCreate(input_grid, &
                            typekind=ESMF_TYPEKIND_R8, &
                            meshloc=ESMF_MESHLOC_ELEMENT, &
                            name=input_diag_names(i), rc=rc)
          if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
          call error_handler("IN FieldCreate", rc)
        endif
    enddo

    input_diag_bundle = ESMF_FieldBundleCreate(fieldList=diag_fields, &
                                    name="input diag data", rc=rc)
    if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldBundleCreate", rc)

    deallocate(diag_fields, input_diag_names)
 end subroutine init_input_diag_fields

!> Read input grid hist data.
!!
!! @param[in] localpet  ESMF local persistent execution thread
!! @author Larissa Reames CIWRO/NOAA/NSSL
 subroutine read_input_hist_data(localpet, npets, mapping)

 character(len=500)              :: the_file
 character(len=50)               :: vname, att_text

 integer, intent(in)             :: localpet, npets, mapping(:)
 integer                         :: error, ncid, rc
 integer                         :: id_dim
 integer                         :: id_var, i, j, nodes

 type(esmf_field),allocatable, target    :: fields(:)

 real(esmf_kind_r8), allocatable, target :: dummy1(:), dummy2(:,:), dummy3(:,:,:), scatter1(:), scatter2(:,:)

 real(esmf_kind_r8), pointer     :: varptr(:), varptr2(:,:)


 the_file = trim(hist_file_input_grid)
 error=nf90_open(trim(the_file),nf90_nowrite,ncid)
 call netcdf_err(error, 'opening: '//trim(the_file) )

 !nSoilLevels must be read before init_input_hist_fields
 if (localpet==0) print*,'- READ nSoilLevels'
 error = nf90_inq_dimid(ncid,'nSoilLevels',id_dim)
 call netcdf_err(error, 'reading nSoilLevels id')

 error=nf90_inquire_dimension(ncid,id_dim,len=nsoil_input)
 call netcdf_err(error, 'reading nSoilLevels')

 allocate(zs_target_grid(nsoil_input,1))
 allocate(dzs_target_grid(nsoil_input,1))

 call init_input_hist_fields(localpet)

 if (localpet==0) print*,"- READ INPUT HIST DATA."

!---------------------------------------------------------------------------
! Read global attributes for use when creating output file
!---------------------------------------------------------------------------

 if (localpet==0) print*,'- READ GLOBAL ATTRIBUTE LSM SCHEME'
 error = nf90_get_att(ncid,NF90_GLOBAL,'config_lsm_scheme',att_text)
 !call netcdf_err(error, 'reading config_lsm_scheme')
 if (error .ne. 0) then
    lsm_scheme = 0
 elseif (INDEX(att_text,'noah') .ne. 0)  then
    lsm_scheme = 2
 elseif (INDEX(att_text, 'ruc') .ne. 0)  then
    lsm_scheme = 3
 endif

 if (localpet==0) print*,'- READ GLOBAL ATTRIBUTE START TIME'
 error = nf90_get_att(ncid,NF90_GLOBAL,'config_start_time',start_time)
 call netcdf_err(error, 'reading config_start_time')

 if (localpet==0) print*,'- READ GLOBAL ATTRIBUTE LMP SCHEME'
 error = nf90_get_att(ncid,NF90_GLOBAL,'config_microp_scheme',att_text)
 !call netcdf_err(error, 'reading config_microp_scheme')
  if (error .ne. 0) then
     mp_scheme = 0
  elseif (trim(att_text) == 'mp_thompson') then
    mp_scheme = 8
  elseif (trim(att_text) == 'mp_nssl2m') then
    mp_scheme = 18
  endif


 if (localpet==0) print*,'- READ GLOBAL ATTRIBUTE CONVECTION SCHEME'
 error = nf90_get_att(ncid,NF90_GLOBAL,'config_convection_scheme',att_text)
 !call netcdf_err(error, 'reading config_conv_scheme')
  if (error .ne. 0) then
    conv_scheme = 0
  elseif (trim(att_text) == 'cu_ntiedke') then
    conv_scheme = 16
  elseif (trim(att_text) == 'cu_kain_fritsch') then
    conv_scheme = 1
  elseif (trim(att_text) == 'cu_grell_freitas') then
    conv_scheme = 3
  endif

 if (localpet==0) print*,'- READ GLOBAL ATTRIBUTE CONFIG_DT'
 error = nf90_get_att(ncid,NF90_GLOBAL,'config_dt',config_dt)
 !call netcdf_err(error, 'reading config_dt')
 if (error .ne. 0) config_dt = 0.0

 vname = 'xtime'
 if (localpet==0) print*, '- READ TIMES VARIABLE'
 error = nf90_inq_dimid(ncid,'StrLen',id_var)
 call netcdf_err(error, 'reading strlen dim id')
 error = nf90_inquire_dimension(ncid,id_var,len=strlen)
 if ( .not. allocated(valid_time)) allocate(valid_time(1,strlen)) !may have been allocated in read_input_diag_data, hence the condition

 error = nf90_inq_varid(ncid,vname,id_var)
 call netcdf_err(error, 'reading xtime id')
 if (localpet==0) print*, "- GETTING xtime"
 error = nf90_get_var(ncid, id_var, valid_time)
 call netcdf_err(error, 'getting xtime')

 ! SOIL CENTER DEPTHS
 if (localpet==0) print*,'- READ ZS ID'
 error=nf90_inq_varid(ncid, 'zs', id_var)
 call netcdf_err(error, 'reading zs id')

 if (localpet==0) print*,'- READ ZS'
 error=nf90_get_var(ncid, id_var, start=(/1,1/),count=(/nsoil_input,1/),values=zs_target_grid)
 call netcdf_err(error, 'reading ZS')

 ! SOIL LAYER THICKNESSES
 if (localpet==0) print*,'- READ DZS ID'
 error=nf90_inq_varid(ncid, 'dzs', id_var)
 call netcdf_err(error, 'reading dzs id')

 if (localpet==0) print*,'- READ DZS'
 error=nf90_get_var(ncid, id_var, start=(/1,1/),count=(/nsoil_input,1/),values=dzs_target_grid)
 call netcdf_err(error, 'reading DZS')

!---------------------------------------------------------------------------
! Initialize 2d esmf atmospheric fields for bilinear/patch interpolation
!---------------------------------------------------------------------------

 if (localpet==0) print*, "Begin reading variables"
 if (n_hist_fields_2d_patch > 0) then
    if (localpet==0) print*, "Read 2d history variables"
    allocate(fields(n_hist_fields_2d_patch))
    allocate(target_hist_units_2d_patch(n_hist_fields_2d_patch))
    allocate(target_hist_longname_2d_patch(n_hist_fields_2d_patch))
     call ESMF_FieldBundleGet(input_hist_bundle_2d_patch, fieldList=fields, &
                              itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                              rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldBundleGet", rc)

     call ESMF_MeshGet(input_grid, nodeCount = nodes, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN MeshGet", rc)

     if(localpet == 0) then
        allocate(dummy1(nCells_input))
        allocate(scatter1(nCells_input))
     else
        allocate(dummy1(1))
        allocate(scatter1(1))
     endif

     do i = 1,n_hist_fields_2d_patch

        call ESMF_FieldGet(fields(i), name=vname, rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)
        call ESMF_FieldGet(fields(i), farrayPtr=varptr, rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

        if (localpet==0) print*,"- READ ", trim(vname)
        error=nf90_inq_varid(ncid, trim(vname), id_var)
        call netcdf_err(error, 'reading field id' )

        call scatter_ncvar_2d(fields(i), localpet, dummy1, scatter1, mapping, ncid, id_var)

        error=nf90_get_att(ncid,id_var,'units',target_hist_units_2d_patch(i))
        call netcdf_err(error, 'reading field units' )
        error=nf90_get_att(ncid,id_var,'long_name',target_hist_longname_2d_patch(i))
        call netcdf_err(error, 'reading field long_name' )

        nullify(varptr)
     enddo
     deallocate(dummy1)
     deallocate(scatter1)
     deallocate(fields)
 endif

!---------------------------------------------------------------------------
! Initialize 2d esmf atmospheric fields for conservative interpolation
!---------------------------------------------------------------------------

 if (n_hist_fields_2d_cons > 0) then
    if (localpet==0) print*, "Read 2d history cons variables"
    allocate(fields(n_hist_fields_2d_cons))
    allocate(target_hist_units_2d_cons(n_hist_fields_2d_cons))
    allocate(target_hist_longname_2d_cons(n_hist_fields_2d_cons))
     call ESMF_FieldBundleGet(input_hist_bundle_2d_cons, fieldList=fields, &
                              itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                              rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldBundleGet", rc)

     call ESMF_MeshGet(input_grid, nodeCount = nodes, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN MeshGet", rc)

     if(localpet == 0) then
        allocate(dummy1(nCells_input))
        allocate(scatter1(nCells_input))
     else
        allocate(dummy1(1))
        allocate(scatter1(1))
     endif

     do i = 1,n_hist_fields_2d_cons

        call ESMF_FieldGet(fields(i), name=vname, rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

        call ESMF_FieldGet(fields(i), farrayPtr=varptr, rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

        if (localpet==0) print*,"- READ ", trim(vname)
        error=nf90_inq_varid(ncid, trim(vname), id_var)
        call netcdf_err(error, 'reading field id' )

        call scatter_ncvar_2d(fields(i), localpet, dummy1, scatter1, mapping, ncid, id_var)

        error=nf90_get_att(ncid,id_var,'units',target_hist_units_2d_cons(i))
        call netcdf_err(error, 'reading field units' )
        error=nf90_get_att(ncid,id_var,'long_name',target_hist_longname_2d_cons(i))
        call netcdf_err(error, 'reading field long_name' )
     enddo
     deallocate(dummy1)
     deallocate(scatter1)
     deallocate(fields)
 endif

!---------------------------------------------------------------------------------
! Initialize 2d esmf atmospheric fields for nearest source to dest interpolation
!---------------------------------------------------------------------------------

 if (n_hist_fields_2d_nstd > 0) then
    if (localpet==0) print*, "Read 2d history nstd variables"
    allocate(fields(n_hist_fields_2d_nstd))
    allocate(target_hist_units_2d_nstd(n_hist_fields_2d_nstd))
    allocate(target_hist_longname_2d_nstd(n_hist_fields_2d_nstd))
     call ESMF_FieldBundleGet(input_hist_bundle_2d_nstd, fieldList=fields, &
                              itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                              rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldBundleGet", rc)

     call ESMF_MeshGet(input_grid, nodeCount = nodes, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN MeshGet", rc)

     if(localpet == 0) then
        allocate(dummy1(nCells_input))
        allocate(scatter1(nCells_input))
     else
        allocate(dummy1(1))
        allocate(scatter1(1))
     endif

     do i = 1,n_hist_fields_2d_nstd

        call ESMF_FieldGet(fields(i), name=vname, rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

        call ESMF_FieldGet(fields(i), farrayPtr=varptr, rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

        if (localpet==0) print*,"- READ ", trim(vname)
        error=nf90_inq_varid(ncid, trim(vname), id_var)
        call netcdf_err(error, 'reading field id' )

        call scatter_ncvar_2d(fields(i), localpet, dummy1, scatter1, mapping, ncid, id_var)

        error=nf90_get_att(ncid,id_var,'units',target_hist_units_2d_nstd(i))
        call netcdf_err(error, 'reading field units' )
        error=nf90_get_att(ncid,id_var,'long_name',target_hist_longname_2d_nstd(i))
        call netcdf_err(error, 'reading field long_name' )
     enddo
     deallocate(dummy1)
     deallocate(scatter1)
     deallocate(fields)
 endif

!---------------------------------------------------------------------------
! Initialize 3d esmf soil fields for bilinear/patch interpolation
!---------------------------------------------------------------------------

 if (n_hist_fields_soil > 0) then
    if (localpet==0) print*, "Read 3d history soil variables"
    allocate(fields(n_hist_fields_soil))
    allocate(target_hist_units_soil(n_hist_fields_soil))
    allocate(target_hist_longname_soil(n_hist_fields_soil))
     call ESMF_FieldBundleGet(input_hist_bundle_soil, fieldList=fields, &
                              itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                              rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldBundleGet", rc)

     call ESMF_MeshGet(input_grid, nodeCount = nodes, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN MeshGet", rc)

     if(localpet == 0) then
        allocate(dummy3(nsoil_input,nCells_input,1))
        allocate(scatter2(nCells_input,nsoil_input))
     else
        allocate(dummy3(1,1,1))
        allocate(scatter2(1,1))
     endif

     do i = 1,n_hist_fields_soil

        call ESMF_FieldGet(fields(i), name=vname, rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

        call ESMF_FieldGet(fields(i), farrayPtr=varptr2, rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

        if (localpet==0) print*,"- READ ", trim(vname)
        error=nf90_inq_varid(ncid, trim(vname), id_var)
        call netcdf_err(error, 'reading field id' )

        call scatter_ncvar_3d(fields(i), localpet, dummy3, scatter2, mapping, ncid, id_var, nsoil_input)

        error=nf90_get_att(ncid,id_var,'units',target_hist_units_soil(i))
        call netcdf_err(error, 'reading field units' )
        error=nf90_get_att(ncid,id_var,'long_name',target_hist_longname_soil(i))
        call netcdf_err(error, 'reading field long_name' )
     enddo
     deallocate(dummy3)
     deallocate(scatter2)
     deallocate(fields)
 endif


!---------------------------------------------------------------------------
! Initialize 3d esmf atmospheric fields with nVertLevels vertical dimension
!---------------------------------------------------------------------------

 if (n_hist_fields_3d_nz > 0 ) then
    if (localpet==0) print*, "Read 3d history atmosphere nz variables"
     allocate(fields(n_hist_fields_3d_nz))
     allocate(target_hist_units_3d_nz(n_hist_fields_3d_nz))
    allocate(target_hist_longname_3d_nz(n_hist_fields_3d_nz))
     call ESMF_FieldBundleGet(input_hist_bundle_3d_nz, fieldList=fields, &
                              itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                              rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldBundleGet", rc)

     call ESMF_MeshGet(input_grid, nodeCount = nodes, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN MeshGet", rc)

     if(localpet == 0) then
        allocate(dummy3(nz_input,nCells_input,1))
        allocate(scatter2(nCells_input,nz_input))
     else
        allocate(dummy3(1,1,1))
        allocate(scatter2(1,1))
     endif

     do i = 1,n_hist_fields_3d_nz

        call ESMF_FieldGet(fields(i), name=vname, rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

        call ESMF_FieldGet(fields(i), farrayPtr=varptr2, rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

        if (localpet==0) print*,"- READ ", trim(vname)
        error=nf90_inq_varid(ncid, trim(vname), id_var)
        call netcdf_err(error, 'reading field id' )

        call scatter_ncvar_3d(fields(i), localpet, dummy3, scatter2, mapping, ncid, id_var, nz_input)

        error=nf90_get_att(ncid,id_var,'units',target_hist_units_3d_nz(i))
        call netcdf_err(error, 'reading field units' )
        error=nf90_get_att(ncid,id_var,'long_name',target_hist_longname_3d_nz(i))
        call netcdf_err(error, 'reading field long_name' )
     enddo
     deallocate(dummy3)
     deallocate(scatter2)
     deallocate(fields)
 endif

!---------------------------------------------------------------------------
! Initialize 3d esmf atmospheric fields U and V if requested
!---------------------------------------------------------------------------
 if (do_u_interp==1) then
    if(localpet == 0) then
       allocate(dummy3(nz_input,nCells_input,1))
       allocate(scatter2(nCells_input,nz_input))
    else
       allocate(dummy3(1,1,1))
       allocate(scatter2(1,1))
    endif

    call ESMF_fieldGet(u_input_grid, farrayPtr=varptr2,rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)
    if (localpet==0) print*, "- READ uReconstructZonal "
    error=nf90_inq_varid(ncid, "uReconstructZonal", id_var)
    call netcdf_err(error, 'reading field id' )

    call scatter_ncvar_3d(u_input_grid, localpet, dummy3, scatter2, mapping, ncid, id_var, nz_input)
    deallocate(dummy3)
    deallocate(scatter2)
 endif

  if (do_v_interp==1) then
     if(localpet == 0) then
        allocate(dummy3(nz_input,nCells_input,1))
        allocate(scatter2(nCells_input,nz_input))
     else
        allocate(dummy3(1,1,1))
        allocate(scatter2(1,1))
     endif

     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)
    if (localpet==0) print*, "- READ uReconstructMeridional "
    error=nf90_inq_varid(ncid, "uReconstructMeridional", id_var)
    call netcdf_err(error, 'reading field id' )
    call scatter_ncvar_3d(v_input_grid, localpet, dummy3, scatter2, mapping, ncid, id_var, nz_input)
    deallocate(dummy3)
    deallocate(scatter2)
 endif

!-------------------------------------------------------------------------------
! Initialize 3d esmf atmospheric fields with nVertLevels+1 vertical dimension
!-------------------------------------------------------------------------------

 if (n_hist_fields_3d_nzp1 > 0 ) then
    if (localpet==0) print*, "Read 3d history atmosphere nz+1 variables"
     allocate(fields(n_hist_fields_3d_nzp1))
     allocate(target_hist_units_3d_nzp1(n_hist_fields_3d_nzp1))
     allocate(target_hist_longname_3d_nzp1(n_hist_fields_3d_nzp1))
     call ESMF_FieldBundleGet(input_hist_bundle_3d_nzp1, fieldList=fields, &
                              itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                              rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldBundleGet", rc)

     call ESMF_MeshGet(input_grid, nodeCount = nodes, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN MeshGet", rc)

     if(localpet == 0) then
        allocate(dummy3(nzp1_input,nCells_input,1))
        allocate(scatter2(nCells_input,nzp1_input))
     else
        allocate(dummy3(1,1,1))
        allocate(scatter2(1,1))
     endif

     do i = 1,n_hist_fields_3d_nzp1

        call ESMF_FieldGet(fields(i), name=vname, rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

        call ESMF_FieldGet(fields(i), farrayPtr=varptr2, rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldGet", rc)

        if (localpet==0) print*,"- READ ", trim(vname)
        error=nf90_inq_varid(ncid, trim(vname), id_var)
        call netcdf_err(error, 'reading field id' )
        error=nf90_get_var(ncid, id_var, dummy3)
        call netcdf_err(error, 'reading field' )

        call scatter_ncvar_3d(fields(i),localpet, dummy3, scatter2, mapping, ncid, id_var, nzp1_input)

        error=nf90_get_att(ncid,id_var,'long_name',target_hist_longname_3d_nzp1(i))
        call netcdf_err(error, 'reading field long_name' )
     enddo
     deallocate(dummy3)
     deallocate(scatter2)
     deallocate(fields)
 endif

!---------------------------------------------------------------------------
! Initialize 3d esmf atmospheric fields with nVertLevels vertical dimension and
! defined on unstructured mesh nodes
!---------------------------------------------------------------------------

 if (n_hist_fields_3d_vert > 0 ) then
    if (localpet==0) print*, "Read 3d history unstructured mesh variables"
     allocate(fields(n_hist_fields_3d_vert))
     allocate(target_hist_units_3d_vert(n_hist_fields_3d_vert))
    allocate(target_hist_longname_3d_vert(n_hist_fields_3d_vert))
     call ESMF_FieldBundleGet(input_hist_bundle_3d_vert, fieldList=fields, &
                              itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                              rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
         call error_handler("IN FieldBundleGet", rc)

     call ESMF_MeshGet(input_grid, nodeCount = nodes, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
         call error_handler("IN MeshGet", rc)

     if(localpet == 0) then
        allocate(dummy3(nz_input,nVert_input,1))
        allocate(scatter2(nVert_input,nz_input))
     else
        allocate(dummy3(1,1,1))
        allocate(scatter2(1,1))
     endif

     do i = 1,n_hist_fields_3d_vert

        call ESMF_FieldGet(fields(i), farrayPtr=varptr2, rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
         call error_handler("IN FieldGet", rc)

        if (localpet==0) print*,"- READ ", trim(vname)
        error=nf90_inq_varid(ncid, trim(vname), id_var)
        call netcdf_err(error, 'reading field id' )

        call scatter_ncvar_3d(fields(i), localpet, dummy3, scatter2, mapping, ncid, id_var, nz_input)

        error=nf90_get_att(ncid,id_var,'units',target_hist_units_3d_vert(i))
        call netcdf_err(error, 'reading field units' )
        error=nf90_get_att(ncid,id_var,'long_name',target_hist_longname_3d_vert(i))
        call netcdf_err(error, 'reading field long_name' )
     enddo
     deallocate(dummy3)
     deallocate(scatter2)
     deallocate(fields)
 endif


 error = nf90_close(ncid)


 end subroutine read_input_hist_data

 subroutine init_input_hist_fields(localpet)

    implicit none

    integer, intent(in)           :: localpet
    integer                       :: i, j, k, n, rc
    integer                       :: n_hist_fields_2d, n_hist_fields_3d
    !type(esmf_field), allocatable :: hist_fields_2d_cons(:), &
 !                                 hist_fields_2d_nstd(:), &
 !                                 hist_fields_2d_patch(:)
    character(50), allocatable    :: input_hist_names_2d(:), &
                                     input_hist_names_2d_cons(:), &
                                     input_hist_names_2d_nstd(:), &
                                     input_hist_names_2d_patch(:), &
                                     input_hist_names_3d(:), &
                                     input_hist_names_3d_nz(:), &
                                     input_hist_names_soil(:), &
                                     input_hist_names_3d_nzp1(:), &
                                     input_hist_names_3d_vert(:), &
                                     target_hist_names_2d(:), &
                                     target_hist_names_3d(:)
    type(esmf_field), allocatable :: fields(:)
    character(50)                 :: cons_vars(2), nstd_vars(4), nzp1_vars(2), &
                                     vert_vars(1)
    character(50)                 :: fname

    cons_vars = (/'snow ','snowh'/)
    nstd_vars = (/'ivgtyp  ','isltyp  ', 'xland   ','landmask'/)
    nzp1_vars = (/'zgrid','w    '/)
    vert_vars = (/'vorticity'/)
    n_hist_fields_2d_cons = 0
    n_hist_fields_2d_nstd = 0
    n_hist_fields_2d_patch = 0
    n_hist_fields_3d_vert = 0
    n_hist_fields_3d_nz = 0
    n_hist_fields_3d_nzp1 = 0

    fname = 'histlist_2d'
    call read_varlist(localpet, fname,n_hist_fields_2d,input_hist_names_2d, target_hist_names_2d)
    fname = 'histlist_3d'
    call read_varlist(localpet,fname,n_hist_fields_3d,input_hist_names_3d, target_hist_names_3d)
    fname = 'histlist_soil'
    call read_varlist(localpet,fname,n_hist_fields_soil,input_hist_names_soil, target_hist_names_soil)

    do i = 1, n_hist_fields_2d
        if (any(cons_vars == input_hist_names_2d(i))) then
            n_hist_fields_2d_cons = n_hist_fields_2d_cons + 1
        elseif (any(nstd_vars == input_hist_names_2d(i))) then
            n_hist_fields_2d_nstd = n_hist_fields_2d_nstd + 1
        else
            n_hist_fields_2d_patch = n_hist_fields_2d_patch + 1
        endif
    enddo

    allocate(input_hist_names_2d_cons(n_hist_fields_2d_cons))
    allocate(input_hist_names_2d_nstd(n_hist_fields_2d_nstd))
    allocate(input_hist_names_2d_patch(n_hist_fields_2d_patch))
    allocate(target_hist_names_2d_cons(n_hist_fields_2d_cons))
    allocate(target_hist_names_2d_nstd(n_hist_fields_2d_nstd))
    allocate(target_hist_names_2d_patch(n_hist_fields_2d_patch))


    j = 0
    k = 0
    n = 0
    do i = 1, n_hist_fields_2d
        if (any(cons_vars == input_hist_names_2d(i))) then
            j = j+1
            input_hist_names_2d_cons(j) = input_hist_names_2d(i)
            target_hist_names_2d_cons(j) = target_hist_names_2d(i)

        elseif (any(nstd_vars == input_hist_names_2d(i))) then
            k = k+1
            input_hist_names_2d_nstd(k) = input_hist_names_2d(i)
            target_hist_names_2d_nstd(k) = target_hist_names_2d(i)
        else
            n = n+1
            input_hist_names_2d_patch(n) = input_hist_names_2d(i)
            target_hist_names_2d_patch(n) = target_hist_names_2d(i)
        endif
    enddo

    do i = 1, n_hist_fields_3d
        if(localpet==0) print*, i
        if (wrf_mod_vars .and. input_hist_names_3d(i)=="uReconstructZonal") then
           do_u_interp = 1
           cycle
        elseif ( wrf_mod_vars .and. input_hist_names_3d(i)=="uReconstructMeridional") then
           do_v_interp = 1
           cycle
        elseif (any(nzp1_vars == input_hist_names_3d(i))) then
            n_hist_fields_3d_nzp1 = n_hist_fields_3d_nzp1 + 1
        elseif (any(vert_vars == input_hist_names_3d(i))) then
            n_hist_fields_3d_vert = n_hist_fields_3d_vert + 1
        else
            n_hist_fields_3d_nz = n_hist_fields_3d_nz + 1
        endif
    enddo

    allocate(input_hist_names_3d_nz(n_hist_fields_3d_nz))
    allocate(input_hist_names_3d_nzp1(n_hist_fields_3d_nzp1))
    allocate(input_hist_names_3d_vert(n_hist_fields_3d_vert))
    allocate(target_hist_names_3d_nz(n_hist_fields_3d_nz))
    allocate(target_hist_names_3d_nzp1(n_hist_fields_3d_nzp1))
    allocate(target_hist_names_3d_vert(n_hist_fields_3d_vert))

    j = 0
    k = 0
    n = 0
    if (do_u_interp==1) then
       if (localpet==0) print*, "- INIT FIELD uReconstructZonal"
       u_input_grid = ESMF_FieldCreate(input_grid, &
                            typekind=ESMF_TYPEKIND_R8, &
                            meshloc=ESMF_MESHLOC_ELEMENT, &
                            name="uReconstructZonal", &
                            ungriddedLBound=(/1/), &
                            ungriddedUBound=(/nz_input/), rc=rc)
       if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldCreate", rc)
    endif

    if (do_v_interp==1) then
       if (localpet==0) print*, "- INIT FIELD uReconstructMeridonal"
       v_input_grid = ESMF_FieldCreate(input_grid, &
                            typekind=ESMF_TYPEKIND_R8, &
                            meshloc=ESMF_MESHLOC_ELEMENT, &
                            name="uReconstructMeridonal", &
                            ungriddedLBound=(/1/), &
                            ungriddedUBound=(/nz_input/), rc=rc)
       if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldCreate", rc)
    endif
    do i = 1, n_hist_fields_3d
        if(localpet==0) print*, i, input_hist_names_3d(i)
        if (do_u_interp==1 .and. input_hist_names_3d(i)=="uReconstructZonal") then
           cycle
        elseif (do_v_interp==1 .and. input_hist_names_3d(i) == "uReconstructMeridional") then
           cycle
        elseif (any(nzp1_vars == input_hist_names_3d(i))) then
            j = j+1
            input_hist_names_3d_nzp1(j) = input_hist_names_3d(i)
            target_hist_names_3d_nzp1(j) = target_hist_names_3d(i)
        elseif (any(vert_vars == input_hist_names_3d(i))) then
            n = n + 1
            input_hist_names_3d_vert(n) = input_hist_names_3d(i)
            target_hist_names_3d_vert(n) = target_hist_names_3d(i)
        else
            k = k+1
            input_hist_names_3d_nz(k) = input_hist_names_3d(i)
            target_hist_names_3d_nz(k) = target_hist_names_3d(i)

        endif
    enddo


    if (localpet==0) print*,"- INITIALIZE INPUT HIST FIELDS."
    if (n_hist_fields_2d_cons > 0) then
        allocate(fields(n_hist_fields_2d_cons))
        do i = 1, n_hist_fields_2d_cons

            if (localpet==0) print*, "- INIT 2D CONSERVATIVE FIELD ", input_hist_names_2d_cons(i)

            fields(i) = ESMF_FieldCreate(input_grid, &
                                typekind=ESMF_TYPEKIND_R8, &
                                meshloc=ESMF_MESHLOC_ELEMENT, &
                                name=input_hist_names_2d_cons(i), rc=rc)
            if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldCreate", rc)
        enddo


        input_hist_bundle_2d_cons = ESMF_FieldBundleCreate(fieldList=fields, &
                                        name="input hist 2d data cons", rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldBundleCreate", rc)
        deallocate(fields)
    endif

    if (n_hist_fields_2d_nstd > 0) then
        allocate(fields(n_hist_fields_2d_nstd))
        do i = 1, n_hist_fields_2d_nstd

            if (localpet==0) print*, "- INIT 2D NEAREST FIELD ", input_hist_names_2d_nstd(i)

            fields(i) = ESMF_FieldCreate(input_grid, &
                                typekind=ESMF_TYPEKIND_R8, &
                                meshloc=ESMF_MESHLOC_ELEMENT, &
                                name=input_hist_names_2d_nstd(i), rc=rc)
            if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldCreate", rc)
        enddo


        input_hist_bundle_2d_nstd = ESMF_FieldBundleCreate(fieldList=fields, &
                                        name="input hist 2d data nstd", rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldBundleCreate", rc)

        deallocate(fields)
    endif

    if (n_hist_fields_2d_patch > 0) then
        allocate(fields(n_hist_fields_2d_patch))
        do i = 1, n_hist_fields_2d_patch

            if (localpet==0) print*, "- INIT 2D BILINEAR FIELD ", input_hist_names_2d_patch(i)

            fields(i) = ESMF_FieldCreate(input_grid, &
                                typekind=ESMF_TYPEKIND_R8, &
                                meshloc=ESMF_MESHLOC_ELEMENT, &
                                name=input_hist_names_2d_patch(i), rc=rc)
            if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldCreate", rc)
        enddo


        input_hist_bundle_2d_patch = ESMF_FieldBundleCreate(fieldList=fields, &
                                        name="input hist 2d data patch", rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldBundleCreate", rc)

        deallocate(fields)
    endif

    if (n_hist_fields_soil > 0) then
        allocate(fields(n_hist_fields_soil))
        do i = 1, n_hist_fields_soil

            if (localpet==0) print*, "- INIT SOIL FIELD ", input_hist_names_soil(i)

            fields(i) = ESMF_FieldCreate(input_grid, &
                                typekind=ESMF_TYPEKIND_R8, &
                                meshloc=ESMF_MESHLOC_ELEMENT, &
                                name=input_hist_names_soil(i), &
                                ungriddedLBound=(/1/), &
                                ungriddedUBound=(/nsoil_input/), rc=rc)
            if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldCreate", rc)
        enddo

        input_hist_bundle_soil = ESMF_FieldBundleCreate(fieldList=fields, &
                                        name="input hist soil data", rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldBundleCreate", rc)

        deallocate(fields)
    endif

    if (n_hist_fields_3d_nz > 0) then
        allocate(fields(n_hist_fields_3d_nz))
        do i = 1, n_hist_fields_3d_nz

            if (localpet==0) print*, "- INIT 3D NZ FIELD ",i, input_hist_names_3d_nz(i)

            fields(i) = ESMF_FieldCreate(input_grid, &
                                typekind=ESMF_TYPEKIND_R8, &
                                meshloc=ESMF_MESHLOC_ELEMENT, &
                                name=input_hist_names_3d_nz(i), &
                                ungriddedLBound=(/1/), &
                                ungriddedUBound=(/nz_input/), rc=rc)
            if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldCreate", rc)
        enddo

        input_hist_bundle_3d_nz = ESMF_FieldBundleCreate(fieldList=fields, &
                                        name="input hist 3d nz data", rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldBundleCreate", rc)

        deallocate(fields)
    endif

    if (n_hist_fields_3d_nzp1 > 0) then
        allocate(fields(n_hist_fields_3d_nzp1))
        do i = 1, n_hist_fields_3d_nzp1

            if (localpet==0) print*, "- INIT 3D NZ+1 FIELD ",i, input_hist_names_3d_nzp1(i)

            fields(i) = ESMF_FieldCreate(input_grid, &
                                typekind=ESMF_TYPEKIND_R8, &
                                meshloc=ESMF_MESHLOC_ELEMENT, &
                                name=input_hist_names_3d_nzp1(i), &
                                ungriddedLBound=(/1/), &
                                ungriddedUBound=(/nzp1_input/), rc=rc)
            if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldCreate", rc)
        enddo

        input_hist_bundle_3d_nzp1 = ESMF_FieldBundleCreate(fieldList=fields, &
                                        name="input hist 3d nzp1 data", rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldBundleCreate", rc)

        deallocate(fields)
    endif

    if (n_hist_fields_3d_vert > 0) then
        allocate(fields(n_hist_fields_3d_vert))
        do i = 1, n_hist_fields_3d_vert

            if (localpet==0) print*, "- INIT 3D VERTICES FIELD ",i,  input_hist_names_3d_vert(i)

            fields(i) = ESMF_FieldCreate(input_grid, &
                                typekind=ESMF_TYPEKIND_R8, &
                                meshloc=ESMF_MESHLOC_NODE, &
                                name=input_hist_names_3d_vert(i), &
                                ungriddedLBound=(/1/), &
                                ungriddedUBound=(/nz_input/), rc=rc)
            if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldCreate", rc)
        enddo

        input_hist_bundle_3d_vert = ESMF_FieldBundleCreate(fieldList=fields, &
                                        name="input hist 3d vertex data", rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldBundleCreate", rc)

        deallocate(fields)
    endif

    deallocate(input_hist_names_2d, target_hist_names_3d)
    deallocate(input_hist_names_2d_cons)
    deallocate(input_hist_names_2d_nstd)
    deallocate(input_hist_names_2d_patch)
    deallocate(input_hist_names_3d_nz)
    deallocate(input_hist_names_3d_nzp1)
    deallocate(input_hist_names_3d_vert)
    deallocate(input_hist_names_soil)


 end subroutine init_input_hist_fields

subroutine read_varlist(localpet, file,nfields,field_names,field_names_target)

    implicit none

   integer, INTENT(IN)                      :: localpet
   character(50), INTENT(IN)                :: file
   integer, INTENT(OUT)                     :: nfields
   character(50), INTENT(OUT), ALLOCATABLE  :: field_names(:), field_names_target(:)

   integer :: k, istat
   character(200) :: line

    LOGICAL :: file_exist

    INQUIRE( FILE = TRIM(file), EXIST = file_exist )
    IF ( .NOT. file_exist ) THEN
        CALL error_handler("VARLIST FILE "//TRIM(file)//" not exist", 1)
    END IF

   open(14, file=trim(file), form='formatted', iostat=istat)
   if (istat /= 0) then
     call error_handler("OPENING VARLIST FILE", istat)
   endif

   nfields = 0

   !Loop over lines of file to count the number of variables
   do
     read(14, '(A)', iostat=istat) line
     if (istat/=0) exit
     if ( trim(line) .eq. '' ) cycle
     nfields = nfields+1
   enddo

   if (localpet==0) print*, "READING ", nfields, " FIELDS ACCORDING TO ", trim(file)
  !if ( nfields == 0) call error_handler("VARLIST FILE IS EMPTY.", -1) ! ok if there are no fields...

   allocate(field_names(nfields))
   allocate(field_names_target(nfields))

   rewind(14)
    do k = 1,nfields
      read(14, *, iostat=istat) field_names(k), field_names_target(k)
     if (istat /= 0) call error_handler("READING VARLIST FILE", istat)

    enddo
   close(14)

end subroutine read_varlist
 end module input_data
