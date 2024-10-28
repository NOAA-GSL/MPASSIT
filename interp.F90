!> @file
!! @brief Interpolate atmospheric and surface data from MPAS diag, init, or forecast files
!! to a target grid.
!! @author  Larissa Reames CIWRO/NOAA/NSSL

!> Public variables are defined below: "target" indicates field
!! associated with the input grid.
!!
!! @author  Larissa Reames CIWRO/NOAA/NSSL

 module interp

 use esmf
 use netcdf
 use utils_mod
 use misc_definitions_module, only : PROJ_LC
 use program_setup, only          : hist_file_input_grid, &
                                    diag_file_input_grid, &
                                    grid_file_input_grid, &
                                    interp_diag, interp_hist, &
                                    i_target, j_target, &
                                    interp_as_bundle, &
                                    proj_code, stand_lon, &
                                    missing_value

 use model_grid, only             : input_grid, target_grid, &
                                    nCells_input, nVert_input,  &
                                    nz_input, nzp1_input, &
                                    nsoil_input, &
                                    cell_latitude_input_grid, &
                                    cell_longitude_input_grid, &
                                    zgrid_input_grid, &
                                    zgrid_target_grid, &
                                    u_input_grid, &
                                    u_target_grid, &
                                    v_input_grid, &
                                    v_target_grid, &
                                    u_target_grid_nostag, &
                                    v_target_grid_nostag, &
                                    hgt_input_grid, hgt_target_grid, &
                                    target_diag_bundle, &
                                    target_diag_names, &
                                    input_diag_bundle, &
                                    input_hist_bundle_2d_cons, &
                                    input_hist_bundle_2d_nstd, &
                                    input_hist_bundle_2d_patch, &
                                    input_hist_bundle_3d_nz, &
                                    input_hist_bundle_3d_nzp1, &
                                    input_hist_bundle_3d_vert, &
                                    input_hist_bundle_soil, &
                                    target_hist_bundle_2d_patch, &
                                    target_hist_bundle_2d_cons, &
                                    target_hist_bundle_2d_nstd, &
                                    target_hist_bundle_3d_nz, &
                                    target_hist_bundle_3d_nzp1, &
                                    target_hist_bundle_3d_vert, &
                                    target_hist_bundle_soil, &
                                    target_diag_names, &
                                    target_hist_names_2d_cons, &
                                    target_hist_names_2d_nstd, &
                                    target_hist_names_2d_patch, &
                                    target_hist_names_3d_nz, &
                                    target_hist_names_3d_nzp1, &
                                    target_hist_names_3d_vert, &
                                    target_hist_names_soil, &
                                    n_diag_fields, nCellsPerPET, &
                                    n_hist_fields_2d_cons, &
                                    n_hist_fields_2d_nstd, &
                                    n_hist_fields_2d_patch, &
                                    n_hist_fields_3d_nz, &
                                    n_hist_fields_3d_nzp1, &
                                    n_hist_fields_3d_vert, &
                                    n_hist_fields_soil, &
                                    do_u_interp, &
                                    do_v_interp, &
                                    cosa_target_grid, &
                                    sina_target_grid, &
                                    do_u10_interp, &
                                    do_v10_interp, &
                                    u10_ind, v10_ind
 implicit none

 private

 type realptr_2d                         !< array to hold 2d pointers
   real(esmf_kind_r8), pointer :: p(:,:)  !< 2d pointer
 end type realptr_2d

  type realptr_3d                          !< array to hold 3d pointers
   real(esmf_kind_r8), pointer :: p(:,:,:) !< 3d pointer
 end type realptr_3d

 public :: interp_data

 type(esmf_routehandle)           :: rh_patch
 real(esmf_kind_r8), parameter    :: spval = 9.9E10
 integer(esmf_kind_i4), pointer, public   :: unmapped_ptr_bi(:), unmapped_ptr_cons(:), &
                                     unmapped_ptr_nstd(:)
 type(logical)                    :: bilinear_regrid, conservative_regrid, nstd_regrid
 contains

 subroutine interp_data(localpet)

     implicit none

     integer, intent(in)   :: localpet

     if (interp_diag) then
        call interp_diag_data(localpet)
     endif
     if (interp_hist) then
        call interp_hist_data(localpet)
     endif

 end subroutine interp_data

 subroutine fill_missing_bundle(localpet,in_bundle,out_bundle,nd,nx,ny, &
                                method,method_flag,unmapped_ptr)
 
   implicit none
   integer, intent(in)                :: localpet, nd, nx, ny
   type(esmf_fieldbundle),intent(in) :: in_bundle, out_bundle
   type(ESMF_RegridMethod_Flag),intent(in) :: method
   type(logical),intent(inout)        :: method_flag
   integer(esmf_kind_i4), pointer,intent(inout) :: unmapped_ptr(:) 
   type(esmf_field)                   :: field,field2
   type(realptr_2d),allocatable       :: fptr2(:)
   type(realptr_3d),allocatable       :: fptr3(:)
   character(len=50)                  :: fname
   integer                            :: num_fields, i, j, k, ij, l(1), u(1), rc, ndims
   integer                           :: isrctermprocessing
   isrctermprocessing = 1

   if(.not. method_flag) then
      call ESMF_FieldBundleGet(in_bundle,1,field,rc=rc)
      if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
       call error_handler("IN FieldBundleGet", rc)
      call ESMF_FieldBundleGet(out_bundle,1,field2, rc=rc)
      if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
       call error_handler("IN FieldBundleGet", rc)
      call ESMF_FieldRegridStore(field, field2, &
                               regridmethod=method, &
                               routehandle=rh_patch, &
                               srcTermProcessing=isrctermprocessing, &
                               unmappedaction=ESMF_UNMAPPEDACTION_IGNORE,&
                               unmappedDstList=unmapped_ptr,&
                               extrapMethod=ESMF_EXTRAPMETHOD_NONE,&
                               rc=rc)
      method_flag  = .true.
   endif

   l = lbound(unmapped_ptr)
   u = ubound(unmapped_ptr)

   call ESMF_FieldBundleGet(out_bundle,fieldCount=num_fields, rc=rc)
    if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
     call error_handler("IN FieldBundleGet", rc)

   allocate(fptr2(num_fields))
   allocate(fptr3(num_fields))

   do i=1, num_fields
     call ESMF_FieldBundleGet(out_bundle,i,field,rc=rc)
      if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
       call error_handler("IN FieldBundleGet", rc)
     call ESMF_FieldGet(field,name=fname,rc=rc)
       if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
           call error_handler("IN FieldGet", rc)
     if (fname=="REFL_10CM") then
        ndims = 3
     else
        ndims = nd
     endif
     if (ndims==2) then
        if(localpet==0) print*, '- FIELDGET 2D', fname
        call ESMF_FieldGet(field,farrayPtr=fptr2(i)%p,rc=rc)
         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
           call error_handler("IN FieldGet", rc)
     elseif (ndims==3) then
        if(localpet==0) print*, '- FIELDGET 3D', fname
        call ESMF_FieldGet(field,farrayPtr=fptr3(i)%p,rc=rc)
         if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
           call error_handler("IN FieldGet", rc)
     else
        call error_handler("FILL_MISSING ONLY SUPPORTS 2 OR 3 DIMENSIONAL VARIABLES. CHECK ND.", -1)
     endif
   enddo
   do ij = l(1), u(1)
        call ij_to_i_j(unmapped_ptr(ij), nx, ny, i, j)
        do k = 1, num_fields
           call ESMF_FieldBundleGet(out_bundle,k,field,rc=rc)
            if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
              call error_handler("IN FieldBundleGet", rc)
           call ESMF_FieldGet(field,name=fname,rc=rc)
            if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
              call error_handler("IN FieldGet", rc)
           if (fname=="REFL_10CM") then
              ndims = 3
           else
              ndims = nd
           endif
           if(localpet==0) print*, k, ndims
           if (ndims==2) fptr2(k)%p(i,j) = missing_value
           if (ndims==3) fptr3(k)%p(i,j,:) = missing_value 
        enddo
   enddo
   if (allocated(fptr2)) deallocate(fptr2)
   if (allocated(fptr3)) deallocate(fptr3)
 end subroutine fill_missing_bundle

subroutine fill_missing_field(localpet,in_field,out_field,nd,nx,ny,method, &
                              method_flag,unmapped_ptr)

   implicit none
   integer, intent(in)               :: localpet, nd, nx, ny
   type(esmf_field),intent(inout)    :: in_field, out_field
   type(ESMF_RegridMethod_Flag),intent(in) :: method
   type(logical),intent(inout)        :: method_flag
   integer(esmf_kind_i4), pointer,intent(inout) :: unmapped_ptr(:)
   real(esmf_kind_r8),pointer       :: fptr2(:,:)
   real(esmf_kind_r8),pointer       :: fptr3(:,:,:)
   integer                            :: i, j, k, ij, l(1), u(1), rc
   integer                           :: isrctermprocessing
   isrctermprocessing = 1
   if(.not. method_flag) then
      call ESMF_FieldRegridStore(in_field, out_field, &
                               regridmethod=method, &
                               routehandle=rh_patch, &
                               srcTermProcessing=isrctermprocessing, &
                               unmappedaction=ESMF_UNMAPPEDACTION_IGNORE,&
                               unmappedDstList=unmapped_ptr,&
                               rc=rc)
      method_flag  = .true.
   endif

   l = lbound(unmapped_ptr)
   u = ubound(unmapped_ptr)
   if (nd==2) then
      call ESMF_FieldGet(out_field,farrayPtr=fptr2,rc=rc)
       if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
         call error_handler("IN FieldGet", rc)
   elseif (nd==3) then
      call ESMF_FieldGet(out_field,farrayPtr=fptr3,rc=rc)
       if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
         call error_handler("IN FieldGet", rc)
   else
      call error_handler("FILL_MISSING ONLY SUPPORTS 2 OR 3 DIMENSIONAL VARIABLES. CHECK ND.",-1)
   endif
   do ij = l(1), u(1)
        call ij_to_i_j(unmapped_ptr(ij), nx, ny, i, j)
        if (nd==2) fptr2(i,j) = missing_value
        if (nd==3) fptr3(i,j,:) = missing_value
   enddo

 end subroutine fill_missing_field 

 subroutine interp_diag_data(localpet)

    implicit none

    integer, intent(in)              :: localpet
    integer                          :: rc
    integer                          :: isrctermprocessing
    type(ESMF_RegridMethod_Flag)     :: method
    
    call init_target_diag_fields(localpet)

    isrctermprocessing = 1
    method = ESMF_REGRIDMETHOD_BILINEAR

    if (localpet==0) print*,"- CREATE DIAG BUNDLE REGRID ROUTEHANDLE"

    call ESMF_FieldBundleRegridStore(input_diag_bundle, target_diag_bundle, &
                                     regridmethod=method, &
                                     routehandle=rh_patch, &
                                     srcTermProcessing=isrctermprocessing, &
                                     unmappedaction=ESMF_UNMAPPEDACTION_IGNORE,&
                                     rc=rc)

     if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldBundleRegridStore", rc)

     if (localpet==0) print*,"- REGRID DIAG FIELDS "
     call ESMF_FieldBundleRegrid(input_diag_bundle, target_diag_bundle, rh_patch, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldBundleRegrid", rc)

     if (do_u10_interp==1 .and. do_v10_interp==1 .and. proj_code==PROJ_LC) then
        call rotate_winds_cgrid(localpet,2)
     endif

     call fill_missing_bundle(localpet,input_diag_bundle,target_diag_bundle,2,i_target,j_target, &
                               method,bilinear_regrid,unmapped_ptr_bi)

  end subroutine interp_diag_data

  subroutine init_target_diag_fields(localpet)

    implicit none

    integer, intent(in)         :: localpet
    integer                     :: i, rc
    type(esmf_field)            :: diag_fields(n_diag_fields)

    if (localpet==0) print*,"- INITIALIZE TARGET DIAG FIELDS."
    do i = 1, n_diag_fields
        if (target_diag_names(i) == "REFL_10CM") then
          if (localpet==0) print*, "- INIT FIELD ", target_diag_names(i)
          diag_fields(i) = ESMF_FieldCreate(target_grid, &
                            typekind=ESMF_TYPEKIND_R8, &
                            staggerloc=ESMF_STAGGERLOC_CENTER, &
                            ungriddedLBound=(/1/), &
                            ungriddedUBound=(/nz_input/), &
                            name=target_diag_names(i), rc=rc)
          if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__))&
          call error_handler("IN FieldCreate", rc)
        else
          if (localpet==0) print*, "- INIT FIELD ", target_diag_names(i)
          diag_fields(i) = ESMF_FieldCreate(target_grid, &
                            typekind=ESMF_TYPEKIND_R8, &
                            staggerloc=ESMF_STAGGERLOC_CENTER, &
                            name=target_diag_names(i), rc=rc)
          if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldCreate", rc)
        endif
    enddo


    target_diag_bundle = ESMF_FieldBundleCreate(fieldList=diag_fields, &
                                    name="target diag data", rc=rc)
    if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldBundleCreate", rc)


 end subroutine init_target_diag_fields

  subroutine interp_hist_data(localpet)

    implicit none

    integer, intent(in)              :: localpet
    integer                          :: l(1), u(1)
    integer                          :: rc, nfields, ij, i, j, n
    integer                          :: isrctermprocessing
    integer                          :: clb_target(2), cub_target(2)
    type(ESMF_RegridMethod_Flag)     :: method
    type(ESMF_RouteHandle)           :: rh_cons, rh_nstd
    type(ESMF_Field), allocatable    :: fields(:)
    type(ESMF_Field), allocatable    :: fields_input_grid(:), fields_target_grid(:)
    real(esmf_kind_r8), pointer      :: field_ptr2(:,:), field_ptr3(:,:,:)
    integer(esmf_kind_i4), pointer   :: unmapped_ptr_u(:), unmapped_ptr_v(:), &
                                        mask_target_ptr(:,:)
    real(esmf_kind_r8), pointer   :: hgt_target_ptr(:,:)
    logical                          :: u_regrid, v_regrid
    call init_target_hist_fields(localpet)

    isrctermprocessing = 1

     if (localpet==0) print*,"- CREATE HGT PATCH REGRID ROUTEHANDLE"
     method = ESMF_REGRIDMETHOD_BILINEAR
     call ESMF_FieldRegridStore(hgt_input_grid, hgt_target_grid, &
                                         regridmethod=method, &
                                         routehandle=rh_patch, &
                                         srcTermProcessing=isrctermprocessing, &
                                         unmappedaction=ESMF_UNMAPPEDACTION_IGNORE,&
                                         unmappedDstList=unmapped_ptr_bi,& 
                                         rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__))&
            call error_handler("IN FieldBundleRegridStore", rc)

     if (localpet==0) print*,"- PATCH REGRID HGT FIELD "
     call ESMF_FieldRegrid(hgt_input_grid, hgt_target_grid, rh_patch, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__))&
        call error_handler("IN FieldBundleRegrid", rc)
    bilinear_regrid = .true.
    call fill_missing_field(localpet,hgt_input_grid, hgt_target_grid,2,i_target, j_target, &
                            method,bilinear_regrid,unmapped_ptr_bi)

   !-----------------------------------------------------------------------
   ! First, set the mask on the target and input grids.
   !-----------------------------------------------------------------------

    if (localpet==0) print*,"- CALL GridAddItem FOR TARGET GRID."
    call ESMF_GridAddItem(target_grid, &
                       itemflag=ESMF_GRIDITEM_MASK, &
                       staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN GridAddItem", rc)

    if (localpet==0) print*,"- CALL GridGetItem FOR TARGET GRID."
    call ESMF_GridGetItem(target_grid, &
                       itemflag=ESMF_GRIDITEM_MASK, &
                       farrayPtr=mask_target_ptr, rc=rc)
    if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN GridGetItem", rc)

    if (localpet==0) print*,"- CALL FieldGet FOR TARGET GRID HGT."
    call ESMF_FieldGet(hgt_target_grid, &
                    computationalLBound=clb_target, &
                    computationalUBound=cub_target, &
                    farrayPtr=hgt_target_ptr, rc=rc)
    if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGet", rc)

    mask_target_ptr = 0
    where (.not. (hgt_target_ptr ==missing_value)) mask_target_ptr = 1 ! outside input data masked

    if (n_hist_fields_2d_patch > 0) then
        method = ESMF_REGRIDMETHOD_BILINEAR
        if (localpet==0) print*,"- CREATE HIST BUNDLE PATCH REGRID ROUTEHANDLE"

        call ESMF_FieldBundleRegridStore(input_hist_bundle_2d_patch, target_hist_bundle_2d_patch, &
                                         regridmethod=method, &
                                         routehandle=rh_patch, &
                                         dstmaskvalues=(/0/), &
                                         srcTermProcessing=isrctermprocessing, &
                                         unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
                                         rc=rc)

         if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldBundleRegridStore", rc)

        if (localpet==0) print*,"- PATCH REGRID INIT FIELDS "
        call ESMF_FieldBundleRegrid(input_hist_bundle_2d_patch, target_hist_bundle_2d_patch, rh_patch, rc=rc)
         if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
           call error_handler("IN FieldBundleRegrid", rc)

        call fill_missing_bundle(localpet,input_hist_bundle_2d_patch, target_hist_bundle_2d_patch,2,i_target,j_target, &
                                 method,bilinear_regrid,unmapped_ptr_bi)

    endif

    if (n_hist_fields_3d_nz>0) then
      method = ESMF_REGRIDMETHOD_BILINEAR
      call ESMF_FieldBundleRegridStore(input_hist_bundle_3d_nz, target_hist_bundle_3d_nz, &
                                         regridmethod=method, &
                                         routehandle=rh_patch, &
                                         dstmaskvalues=(/0/), &
                                         srcTermProcessing=isrctermprocessing, &
                                         unmappedaction=ESMF_UNMAPPEDACTION_IGNORE,&
                                         rc=rc)

        if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__))&
            call error_handler("IN FieldBundleRegridStore", rc)

      call ESMF_FieldBundleRegrid(input_hist_bundle_3d_nz, target_hist_bundle_3d_nz, rh_patch, rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldBundleRegrid", rc)

      call fill_missing_bundle(localpet,input_hist_bundle_3d_nz, target_hist_bundle_3d_nz,3,i_target,j_target, &
                                 method,bilinear_regrid,unmapped_ptr_bi)

    endif

    if (do_u_interp==1) then
       if (localpet==0) print*, "- CREATE REGRID uReconstructZonal ROUTEHANDLE"
       method = ESMF_REGRIDMETHOD_BILINEAR
       call ESMF_FieldRegridStore(u_input_grid,u_target_grid_nostag, &
                                        regridmethod=method, &
                                        routehandle=rh_patch, &
                                         dstmaskvalues=(/0/), &
                                        srcTermProcessing=isrctermprocessing, &
                                        unmappedaction=ESMF_UNMAPPEDACTION_IGNORE,&
                                        rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
          call error_handler("IN FieldRegridStore", rc)

       call ESMF_FieldRegrid(u_input_grid,u_target_grid_nostag, rh_patch, rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldRegrid", rc)
    endif

    if (do_v_interp==1) then
       if (localpet==0) print*, "- CREATE REGRID uReconstructMeridional ROUTEHANDLE"
       method = ESMF_REGRIDMETHOD_BILINEAR
       call ESMF_FieldRegridStore(v_input_grid,v_target_grid_nostag, &
                                        regridmethod=method, &
                                        routehandle=rh_patch, &
                                         dstmaskvalues=(/0/), &
                                        srcTermProcessing=isrctermprocessing, &
                                        unmappedaction=ESMF_UNMAPPEDACTION_IGNORE,&
                                        rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
          call error_handler("IN FieldRegridStore", rc)

       call ESMF_FieldRegrid(v_input_grid,v_target_grid_nostag, rh_patch, rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldRegrid", rc)
    endif

    if (do_u_interp==1 .and. do_v_interp==1 .and. proj_code==PROJ_LC) then
       call rotate_winds_cgrid(localpet,3)
    endif

    if (do_u_interp==1) then
       if (localpet==0) print*, "- CREATE REGRID uReconstructZonal ROUTEHANDLE"
       method = ESMF_REGRIDMETHOD_BILINEAR
       call ESMF_FieldRegridStore(u_target_grid_nostag, u_target_grid, &
                                        regridmethod=method, &
                                        routehandle=rh_patch, &
                                         dstmaskvalues=(/0/), &
                                        srcTermProcessing=isrctermprocessing, &
                                        unmappedaction=ESMF_UNMAPPEDACTION_IGNORE,&
                                        unmappedDstList=unmapped_ptr_u, &
                                        rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
          call error_handler("IN FieldRegridStore", rc)

       call ESMF_FieldRegrid(u_target_grid_nostag,u_target_grid, rh_patch, rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldRegrid", rc)

       call fill_missing_field(localpet,u_input_grid, u_target_grid,3,i_target+1, j_target, &
                            method,u_regrid,unmapped_ptr_u)     
    endif

    if (do_v_interp==1) then
       if (localpet==0) print*, "- CREATE REGRID uReconstructMeridional ROUTEHANDLE"
       method = ESMF_REGRIDMETHOD_BILINEAR
       call ESMF_FieldRegridStore(v_target_grid_nostag, v_target_grid, &
                                        regridmethod=method, &
                                        routehandle=rh_patch, &
                                         dstmaskvalues=(/0/), &
                                        srcTermProcessing=isrctermprocessing, &
                                        unmappedaction=ESMF_UNMAPPEDACTION_IGNORE,&
                                        rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
          call error_handler("IN FieldRegridStore", rc)

       call ESMF_FieldRegrid(v_target_grid_nostag,v_target_grid, rh_patch, rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldRegrid", rc)
 
       call fill_missing_field(localpet,v_input_grid, v_target_grid,3,i_target, j_target+1, &
                            method,v_regrid,unmapped_ptr_v)
    endif


    if (n_hist_fields_3d_nzp1>0) then
        if (localpet==0) print*,"- CREATE HIST BUNDLE PATCH REGRID ROUTEHANDLE"
       method = ESMF_REGRIDMETHOD_BILINEAR
       call ESMF_FieldBundleRegridStore(input_hist_bundle_3d_nzp1, target_hist_bundle_3d_nzp1, &
                                         regridmethod=method, &
                                         routehandle=rh_patch, &
                                         dstmaskvalues=(/0/), &
                                         srcTermProcessing=isrctermprocessing, &
                                         unmappedaction=ESMF_UNMAPPEDACTION_IGNORE,&
                                         rc=rc)

       if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldBundleRegridStore", rc)

       call ESMF_FieldBundleRegrid(input_hist_bundle_3d_nzp1, target_hist_bundle_3d_nzp1, rh_patch, rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldBundleRegrid", rc)

       call fill_missing_bundle(localpet,input_hist_bundle_3d_nzp1, target_hist_bundle_3d_nzp1,3,i_target,j_target, &
                                 method,bilinear_regrid,unmapped_ptr_bi)
    endif


    if (n_hist_fields_3d_vert>0) then
       if (localpet==0) print*,"- CREATE HIST BUNDLE VERT BILINEAR REGRID ROUTEHANDLE"
       method = ESMF_REGRIDMETHOD_BILINEAR
       call ESMF_FieldBundleRegridStore(input_hist_bundle_3d_vert,target_hist_bundle_3d_vert, &
                                         regridmethod=method, &
                                         routehandle=rh_patch, &
                                         dstmaskvalues=(/0/), &
                                         srcTermProcessing=isrctermprocessing, &
                                         unmappedaction=ESMF_UNMAPPEDACTION_IGNORE,&
                                         rc=rc)

       if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
        call error_handler("IN FieldBundleRegridStore", rc)

       call ESMF_FieldBundleRegrid(input_hist_bundle_3d_vert, target_hist_bundle_3d_vert, rh_patch, rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldBundleRegrid", rc)
       call fill_missing_bundle(localpet,input_hist_bundle_3d_vert, target_hist_bundle_3d_vert,3,i_target,j_target, &
                                 method,bilinear_regrid,unmapped_ptr_bi)
    endif

    if (n_hist_fields_2d_cons>0) then
       if (localpet==0) print*,"- CREATE HIST BUNDLE CONSERVATIVE REGRID ROUTEHANDLE"
       method = ESMF_REGRIDMETHOD_CONSERVE
       if ( interp_as_bundle ) then
          call ESMF_FieldBundleRegridStore(input_hist_bundle_2d_cons, target_hist_bundle_2d_cons, &
                                            regridmethod=method, &
                                            routehandle=rh_cons, &
                                         dstmaskvalues=(/0/), &
                                            srcTermProcessing=isrctermprocessing, &
                                            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE,&
                                            rc=rc)

          if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
           call error_handler("IN FieldBundleRegridStore", rc)

          call ESMF_FieldBundleRegrid(input_hist_bundle_2d_cons, target_hist_bundle_2d_cons, rh_cons, rc=rc)
          if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldBundleRegrid", rc)

          call fill_missing_bundle(localpet,input_hist_bundle_2d_cons, target_hist_bundle_2d_cons,2,i_target,j_target, &
                                 method,conservative_regrid,unmapped_ptr_cons)
       else
          allocate(fields_input_grid(n_hist_fields_2d_cons))
          allocate(fields_target_grid(n_hist_fields_2d_cons))
          call ESMF_FieldBundleGet(input_hist_bundle_2d_cons, fieldList=fields_input_grid, &
                                itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                                rc=rc)
          call ESMF_FieldBundleGet(target_hist_bundle_2d_cons, fieldList=fields_target_grid, &
                                itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                                rc=rc)
          call ESMF_FieldRegridStore(fields_input_grid(1), fields_target_grid(1), &
                                      regridmethod=method, &
                                      routehandle=rh_cons, &
                                         dstmaskvalues=(/0/), &
                                      srcTermProcessing=isrctermprocessing, &
                                      unmappedaction=ESMF_UNMAPPEDACTION_IGNORE,&
                                      unmappedDstList=unmapped_ptr_cons, &
                                      rc=rc)
          if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__))&
             call error_handler("IN FieldRegridStore", rc)
          conservative_regrid = .true.
          do i = 1, n_hist_fields_2d_cons
             call ESMF_FieldRegrid(fields_input_grid(i), fields_target_grid(i), rh_cons, rc=rc)
             if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__))&
                call error_handler("IN FieldRegrid", rc)

             call fill_missing_field(localpet,fields_input_grid(i), fields_target_grid(i),2,i_target, j_target, &
                            method,conservative_regrid,unmapped_ptr_cons)

          enddo

          ! update the fields in target_hist_bundle_2d_cons
          call ESMF_FieldBundleAddReplace(target_hist_bundle_2d_cons, fields_target_grid, rc=rc)
          if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__))&
             call error_handler("IN ESMF_FieldBundleAddReplace", rc)

          deallocate(fields_input_grid,fields_target_grid)
       endif
    endif

    if (n_hist_fields_2d_nstd>0) then
       if (localpet==0) print*,"- CREATE HIST BUNDLE NSTD REGRID ROUTEHANDLE"
       method = ESMF_REGRIDMETHOD_NEAREST_STOD
       call ESMF_FieldBundleRegridStore(input_hist_bundle_2d_nstd, target_hist_bundle_2d_nstd, &
                                         regridmethod=method, &
                                         routehandle=rh_nstd, &
                                         dstmaskvalues=(/0/), &
                                         srcTermProcessing=isrctermprocessing, &
                                         unmappedaction=ESMF_UNMAPPEDACTION_IGNORE,&
                                         extrapMethod=ESMF_EXTRAPMETHOD_NONE,&
                                         rc=rc)

       if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldBundleRegridStore", rc)

       if (localpet==0) print*,"- REGRID HIST BUNDLE NSTD"
       call ESMF_FieldBundleRegrid(input_hist_bundle_2d_nstd, target_hist_bundle_2d_nstd, rh_nstd, rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldBundleRegrid", rc)

       call fill_missing_bundle(localpet,input_hist_bundle_2d_nstd, target_hist_bundle_2d_nstd,2,i_target,j_target, &
                                 ESMF_REGRIDMETHOD_BILINEAR,bilinear_regrid,unmapped_ptr_bi)

    endif

    if (n_hist_fields_soil>0) then
        if (localpet==0) print*,"- CREATE HIST BUNDLE SOIL REGRID ROUTEHANDLE"
        call ESMF_FieldBundleRegridStore(input_hist_bundle_soil, target_hist_bundle_soil, &
                                         regridmethod=method, &
                                         routehandle=rh_nstd, &
                                         dstmaskvalues=(/0/), &
                                         srcTermProcessing=isrctermprocessing, &
                                         unmappedaction=ESMF_UNMAPPEDACTION_IGNORE,&
                                         extrapMethod=ESMF_EXTRAPMETHOD_NONE,&
                                         rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
          call error_handler("IN FieldBundleRegridStore", rc)
   
       call ESMF_FieldBundleRegrid(input_hist_bundle_soil, target_hist_bundle_soil, rh_nstd, rc=rc)
         if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldBundleRegrid", rc)
   
       call fill_missing_bundle(localpet,input_hist_bundle_soil, target_hist_bundle_soil,3,i_target,j_target, &
                                ESMF_REGRIDMETHOD_BILINEAR,bilinear_regrid,unmapped_ptr_bi)                             

    endif

    if (localpet==0) print*,"- CALL FieldRegridRelease."
    call ESMF_FieldBundleRegridRelease(routehandle=rh_patch, rc=rc)
     if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldRegridRelease", rc)
    if(n_hist_fields_2d_cons>0) then
       if (localpet==0) print*,"- CALL FieldRegridRelease."
       call ESMF_FieldBundleRegridRelease(routehandle=rh_cons, rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldRegridRelease", rc)
    endif
    if(n_hist_fields_2d_nstd>0) then
       if (localpet==0) print*,"- CALL FieldRegridRelease."
       call ESMF_FieldBundleRegridRelease(routehandle=rh_nstd, rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldRegridRelease", rc)
    endif
  end subroutine interp_hist_data

  subroutine init_target_hist_fields(localpet)

    implicit none

    integer, intent(in)          :: localpet
    integer                      :: i, n, rc
    character(50), allocatable   :: field_names(:)
    type(esmf_field),allocatable :: fields(:)


    if (localpet==0) print*, "- INIT TARGET HGT FIELD "
    hgt_target_grid = ESMF_FieldCreate(target_grid, &
                              typekind=ESMF_TYPEKIND_R8, &
                              staggerloc=ESMF_STAGGERLOC_CENTER, &
                              name='HGT', rc=rc)
    if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldCreate", rc)

    if (do_u_interp==1) then
         if (localpet==0) print*, "- INIT FIELD U"
         u_target_grid = ESMF_FieldCreate(target_grid, &
                             typekind=ESMF_TYPEKIND_R8, &
                             staggerloc=ESMF_STAGGERLOC_EDGE1, &
                             name="U", & 
                             ungriddedLBound=(/1/), &
                             ungriddedUBound=(/nz_input/),rc=rc)
         if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldCreate", rc)

         u_target_grid_nostag = ESMF_FieldCreate(target_grid, &
                             typekind=ESMF_TYPEKIND_R8, &
                             staggerloc=ESMF_STAGGERLOC_CENTER, &
                             name="UMASS", &
                             ungriddedLBound=(/1/), &
                             ungriddedUBound=(/nz_input/),rc=rc)
         if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldCreate", rc)
    endif
    if(do_v_interp==1) then
         if (localpet==0) print*, "- INIT FIELD V"
         v_target_grid = ESMF_FieldCreate(target_grid, &
                            typekind=ESMF_TYPEKIND_R8, &
                            staggerloc=ESMF_STAGGERLOC_EDGE2, &
                            name="V", &
                             ungriddedLBound=(/1/), &
                             ungriddedUBound=(/nz_input/), rc=rc)
         if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldCreate", rc)

         v_target_grid_nostag = ESMF_FieldCreate(target_grid, &
                            typekind=ESMF_TYPEKIND_R8, &
                            staggerloc=ESMF_STAGGERLOC_CENTER, &
                            name="VMASS", &
                             ungriddedLBound=(/1/), &
                             ungriddedUBound=(/nz_input/), rc=rc)
         if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldCreate", rc)
    endif

    if (localpet==0) print*,"- INITIALIZE TARGET 2D CONS HIST FIELDS."

    if (n_hist_fields_2d_cons > 0) then
        allocate(fields(n_hist_fields_2d_cons))
        do i = 1, n_hist_fields_2d_cons

            if (localpet==0) print*, "- INIT FIELD ", target_hist_names_2d_cons(i)
            fields(i) = ESMF_FieldCreate(target_grid, &
                                typekind=ESMF_TYPEKIND_R8, &
                                staggerloc=ESMF_STAGGERLOC_CENTER, &
                                name=target_hist_names_2d_cons(i), rc=rc)
            if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldCreate", rc)
        enddo


        target_hist_bundle_2d_cons = ESMF_FieldBundleCreate(fieldList=fields, &
                                        name="target init 2d cons data", rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldBundleCreate", rc)
        deallocate(fields)
    endif

    if (localpet==0) print*,"- INITIALIZE TARGET 2D NSTD HIST FIELDS."

    if (n_hist_fields_2d_nstd > 0) then
        allocate(fields(n_hist_fields_2d_nstd))
        do i = 1, n_hist_fields_2d_nstd

            if (localpet==0) print*, "- INIT FIELD ", target_hist_names_2d_nstd(i)
            fields(i) = ESMF_FieldCreate(target_grid, &
                                typekind=ESMF_TYPEKIND_R8, &
                                staggerloc=ESMF_STAGGERLOC_CENTER, &
                                name=target_hist_names_2d_nstd(i), rc=rc)
            if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldCreate", rc)
        enddo


        target_hist_bundle_2d_nstd = ESMF_FieldBundleCreate(fieldList=fields, &
                                        name="target init 2d nstd data", rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldBundleCreate", rc)
        deallocate(fields)
    endif

    if (localpet==0) print*,"- INITIALIZE TARGET 2D PATCH HIST FIELDS."

    if (n_hist_fields_2d_patch > 0) then
        allocate(fields(n_hist_fields_2d_patch))
        do i = 1, n_hist_fields_2d_patch

            if (localpet==0) print*, "- INIT FIELD ", target_hist_names_2d_patch(i)
            fields(i) = ESMF_FieldCreate(target_grid, &
                                typekind=ESMF_TYPEKIND_R8, &
                                staggerloc=ESMF_STAGGERLOC_CENTER, &
                                name=target_hist_names_2d_patch(i), rc=rc)
            if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldCreate", rc)
        enddo


        target_hist_bundle_2d_patch = ESMF_FieldBundleCreate(fieldList=fields, &
                                        name="target init 2d patch data", rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldBundleCreate", rc)
        deallocate(fields)
    endif

    if (n_hist_fields_soil>0) then
        allocate(fields(n_hist_fields_soil))
        do i = 1, n_hist_fields_soil
            if (localpet==0) print*, "- INIT FIELD ", target_hist_names_soil(i)
            fields(i) = ESMF_FieldCreate(target_grid, &
                                typekind=ESMF_TYPEKIND_R8, &
                                staggerloc=ESMF_STAGGERLOC_CENTER, &
                                name=target_hist_names_soil(i), &
                                ungriddedLBound=(/1/), &
                                ungriddedUBound=(/nsoil_input/), rc=rc)
            if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldCreate", rc)
        enddo


        target_hist_bundle_soil= ESMF_FieldBundleCreate(fieldList=fields, &
                                        name="target init soil data", rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldBundleCreate", rc)
        deallocate(fields)
    endif

    if (n_hist_fields_3d_nz>0) then
        allocate(fields(n_hist_fields_3d_nz))
        do i = 1, n_hist_fields_3d_nz
            if (localpet==0) print*, "- INIT FIELD ", target_hist_names_3d_nz(i)
            fields(i) = ESMF_FieldCreate(target_grid, &
                                typekind=ESMF_TYPEKIND_R8, &
                                staggerloc=ESMF_STAGGERLOC_CENTER, &
                                name=target_hist_names_3d_nz(i), &
                                ungriddedLBound=(/1/), &
                                ungriddedUBound=(/nz_input/), rc=rc)
            if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldCreate", rc)
        enddo


        target_hist_bundle_3d_nz = ESMF_FieldBundleCreate(fieldList=fields, &
                                        name="target init 3d nz data", rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldBundleCreate", rc)
        deallocate(fields)
    endif

    if (n_hist_fields_3d_nzp1>0) then
        allocate(fields(n_hist_fields_3d_nzp1))
        do i = 1, n_hist_fields_3d_nzp1
            if (localpet==0) print*, "- INIT FIELD ", target_hist_names_3d_nzp1(i)
            fields(i) = ESMF_FieldCreate(target_grid, &
                                typekind=ESMF_TYPEKIND_R8, &
                                staggerloc=ESMF_STAGGERLOC_CENTER, &
                                name=target_hist_names_3d_nzp1(i), &
                                ungriddedLBound=(/1/), &
                                ungriddedUBound=(/nzp1_input/), rc=rc)
            if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldCreate", rc)
        enddo


        target_hist_bundle_3d_nzp1 = ESMF_FieldBundleCreate(fieldList=fields, &
                                        name="target init 3d nz data", rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
            call error_handler("IN FieldBundleCreate", rc)
        deallocate(fields)
    endif

    if (n_hist_fields_3d_vert>0) then
        allocate(fields(n_hist_fields_3d_vert))
        do i = 1, n_hist_fields_3d_vert
            if (localpet==0) print*, "- INIT FIELD ", target_hist_names_3d_vert(i)
            fields(i) = ESMF_FieldCreate(target_grid, &
                                typekind=ESMF_TYPEKIND_R8, &
                                staggerloc=ESMF_STAGGERLOC_CENTER, &
                                name=target_hist_names_3d_vert(i), &
                                ungriddedLBound=(/1/), &
                                ungriddedUBound=(/nz_input/), rc=rc)
            if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
            call error_handler("IN FieldCreate", rc)
        enddo


        target_hist_bundle_3d_vert = ESMF_FieldBundleCreate(fieldList=fields, &
                                        name="target init 3d nz vert data", rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
            call error_handler("IN FieldBundleCreate", rc)
        deallocate(fields)
    endif

    deallocate(target_hist_names_2d_cons, target_hist_names_2d_nstd)
    deallocate(target_hist_names_2d_patch, target_hist_names_3d_nz)
    deallocate(target_hist_names_3d_nzp1, target_hist_names_3d_vert)

 end subroutine init_target_hist_fields

 subroutine rotate_winds_cgrid(localpet,wind_dim)
      use constants_module
      implicit none
      integer, intent(in) :: localpet, wind_dim
      real(esmf_kind_r8), pointer, dimension(:,:,:) :: u_ptr3, v_ptr3
      real(esmf_kind_r8), pointer, dimension(:,:) :: u_ptr2, v_ptr2
      real(esmf_kind_r8), pointer, dimension(:,:)   :: cosa, sina  
      type(esmf_field), allocatable    :: fields(:) 
      integer, dimension(2)            ::  clb, cub        
      double precision :: tana
      integer :: i,j,rc

      if (wind_dim==3) then
         call ESMF_FieldGet(u_target_grid_nostag,farrayPtr = u_ptr3,rc=rc)
           if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
               call error_handler("IN FieldGet", rc)
         call ESMF_FieldGet(v_target_grid_nostag,farrayPtr = v_ptr3,rc=rc)
           if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
               call error_handler("IN FieldGet", rc)
      elseif(wind_dim==2) then
         allocate (fields(n_diag_fields))
         call ESMF_FieldBundleGet(target_diag_bundle, fieldList=fields, &
                                  itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                                  rc=rc)
         if (ESMF_logFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) &
             call error_handler("IN FieldBundleGet", rc)
         call ESMF_FieldGet(fields(u10_ind),farrayPtr = u_ptr2,rc=rc)
           if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
               call error_handler("IN FieldGet", rc)
         call ESMF_FieldGet(fields(v10_ind),farrayPtr = v_ptr2,rc=rc)
           if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
               call error_handler("IN FieldGet", rc)
         deallocate(fields)
      else
         call error_handler("In rotate_winds_cgrid: input wind_dim must equal 2 (for 10-m winds) or 3 (for 3-d winds)", -1)
      endif

      call ESMF_FieldGet(cosa_target_grid,farrayPtr = cosa, &
                         computationalLBound=clb,&
                         computationalUBound=cub, rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
            call error_handler("IN FieldGet", rc)
      call ESMF_FieldGet(sina_target_grid,farrayPtr = sina, &
                         rc=rc)
        if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
            call error_handler("IN FieldGet", rc)   

      do j = clb(2),cub(2)
      do i = clb(1),cub(1)
         tana = sina(i,j)/cosa(i,j)
         if (wind_dim==3) then
            u_ptr3(i,j,:) = (u_ptr3(i,j,:) + v_ptr3(i,j,:) * tana) / (cosa(i,j) + sina(i,j) * tana)
            v_ptr3(i,j,:) = (v_ptr3(i,j,:) - u_ptr3(i,j,:) * sina(i,j)) / cosa(i,j)
         elseif(wind_dim==2) then
            u_ptr2(i,j) = (u_ptr2(i,j) + v_ptr2(i,j) * tana) / (cosa(i,j) + sina(i,j) * tana)
            v_ptr2(i,j) = (v_ptr2(i,j) - u_ptr2(i,j) * sina(i,j)) / cosa(i,j)
         endif
      enddo
      enddo
   end subroutine rotate_winds_cgrid
   
!> Convert 1d index to 2d indices.
!!
!! @param[in] ij  the 1d index
!! @param[in] itile  i-dimension of the tile
!! @param[in] jtile  j-dimension of the tile
!! @param[out] i  the "i" index
!! @param[out] j  the "j" index
!! @author George Gayno NOAA/EMC
 subroutine ij_to_i_j(ij, itile, jtile, i, j)

 implicit none

 integer(esmf_kind_i4), intent(in)  :: ij
 integer              , intent(in)  :: itile, jtile

 integer              , intent(out) :: i, j

 integer                            :: tile_num
 integer                            :: pt_loc_this_tile

 tile_num = ((ij-1) / (itile*jtile)) ! tile number minus 1
 pt_loc_this_tile = ij - (tile_num * itile * jtile)
                                     ! "ij" location of point within tile.

 j = (pt_loc_this_tile - 1) / itile + 1
 i = mod(pt_loc_this_tile, itile)

 if (i==0) i = itile

 return

 end subroutine ij_to_i_j
end module interp


