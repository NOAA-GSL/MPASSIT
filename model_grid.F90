!> @file
!! @brief Specify input and target model grids.
!! @author Larissa Reames CIWRO/NOAA/NSSL/FRDD

!> Sets up the ESMF grid objects for the input data grid and target grid.
!!
!! @author Larissa Reames CIWRO/NOAA/NSSL/FRDD
 module model_grid

 use esmf
 use ESMF_LogPublicMod
 use utils_mod
 use constants_module
 use misc_definitions_module
 use map_utils_mod
 use program_setup

 implicit none

 private

 integer, public                        :: nCells_input
                                           !< cells on input grid
 integer, public                        :: nVert_input
                                           !< vertices on input grid
 integer, public                        :: nz_input
                                           !< number of input grid atm layers
 integer, public                        :: nzp1_input
                                           !< number of input grid atm layer interfaces
 integer, public                        :: nsoil_input
                                           !< number of input soil levels
 !real, public                           :: dx
 !                                          !< grid size (m) of target grid
 character(50), public                  :: start_time
                                           !< simulation start time input data
 character(20), public,allocatable      :: valid_time(:,:)
                                           !< valid forecast time input data
 integer, public                        :: strlen
                                           !< StrLen on input file
 real, public                           :: config_dt
                                           !< timestep (seconds) of input data
 integer, public                        :: lsm_scheme
                                           !< land surface scheme input data
 integer, public                        :: mp_scheme
                                           !< microphysics scheme input data
 integer, public                        :: conv_scheme
                                           !< convection scheme input data
 !real, public                           :: cen_lat
 !                                          !< target grid projection center latitude
 !real, public                           :: cen_lon
 !                                          !< target grid projection center longitude
 !real, public                           :: moad_cen_lat
 !                                          !< target grid moad center latitude
 !real, public                           :: pole_lat
 !                                          !< target grid projection pole latitude
 !real, public                           :: pole_lon
 !                                          !< target grid projection pole longitude
 integer, public                        :: ip1_target
                                           !< ip1_target plus 1
 integer, public                        :: jp1_target
                                           !< jp1_target plus 1

 integer, allocatable, public           :: elemIDs(:)
                                           !< IDs of the elements on present PET
 integer, allocatable, public           :: nodeIDs(:)
                                            !< IDs of the nodes on present PET
 integer, public                        :: nCellsPerPET
                                            !< Number of cells on this PET 
 integer, public                        :: nNodesPerPET
                                           !< Number of nodes on this PET

 type(esmf_mesh),  public               :: input_grid
                                           !< input grid esmf grid object
 type(esmf_grid),  public               :: target_grid
                                           !< target grid esmf grid object.

 type(esmf_field),  public              :: cell_latitude_input_grid
                                           !< latitude of grid center, input grid
 type(esmf_field),  public              :: cell_longitude_input_grid
                                           !< longitude of grid center, input grid

 type(esmf_field),  public              :: node_latitude_input_grid
                                           !< latitude of grid center, input grid
 type(esmf_field),  public              :: node_longitude_input_grid
                                           !< longitude of grid center, input grid

 type(esmf_field),  public              :: zgrid_input_grid
                                           !< esmf field to hold level height on input grid
 type(esmf_field),  public              :: zgrid_target_grid
                                           !< esmf field to hold level height on target grid
 type(esmf_field),  public              :: u_input_grid
                                           !< esmf field to hold u wind on the input grid
 type(esmf_field),  public              :: u_target_grid
                                           !< esnf field to hold u wind on the target grid
 type(esmf_field),  public              :: v_input_grid
                                           !< esmf field to hold v wind on the input grid
 type(esmf_field),  public              :: v_target_grid
                                           !< esmf field to hold v wind on the target grid
 type(esmf_field),  public              :: u_target_grid_nostag
                                           !< esmf field to hold v wind on the target grid mass points
 type(esmf_field),  public              :: v_target_grid_nostag
                                           !< esmf field to hold v wind on the target grid mass points
 type(esmf_field),  public              :: latitude_target_grid
                                           !< latitude of grid center, target grid
 type(esmf_field),  public              :: longitude_target_grid
                                           !< longitude of grid center, target grid
 type(esmf_field),  public              :: latitude_u_target_grid
                                           !< latitude of grid u stagger, target
                                           !grid
 type(esmf_field),  public              :: longitude_u_target_grid
                                           !< longitude of grid u stagger, target
                                           !grid
 type(esmf_field),  public              :: latitude_v_target_grid
                                           !< latitude of grid v stagger, target
                                           !grid
 type(esmf_field),  public              :: longitude_v_target_grid
                                           !< longitude of grid v stagger, target
                                           !grid
 real(esmf_kind_r8), allocatable , public :: zs_target_grid(:,:)
                                          !< soil center depth, target grid
 type(esmf_field), public               :: hgt_input_grid, hgt_target_grid
                                          !< surface elevation, target grid
 type(esmf_field), public               :: mapfac_m_target_grid
                                          !< target grid map factor at grid
                                          ! center 
 type(esmf_field), public               :: mapfac_u_target_grid
                                          !< target grid map factor at u stagger
                                          ! points
 type(esmf_field), public               :: mapfac_v_target_grid
                                          !< target grid map factor at v stagger
                                          ! points
 type(esmf_field), public               :: sina_target_grid
                                         !< Sine of grid rotation angle
                                         !  valid only for grid_type="lambert"
 type(esmf_field), public                :: cosa_target_grid
                                         !< Cosine of grid rotation angle
                                         !  valid only for grid_type="lambert"
 type(esmf_field), public                :: sina_u_target_grid
                                         !< Sine of grid rotation angle for u winds
                                         !  valid only for grid_type="lambert"
 type(esmf_field), public                :: cosa_u_target_grid
                                         !< Cosine of grid rotation angle for u wind
                                         !  valid only for grid_type="lambert" 
 type(esmf_field), public                :: sina_v_target_grid
                                         !< Sine of grid rotation angle for v wind
                                         !  valid only for grid_type="lambert"
 type(esmf_field), public                :: cosa_v_target_grid
                                         !< Cosine of grid rotation angle for v wind
                                         !  valid only for grid_type="lambert"


 integer, public                       :: n_diag_fields
                                          !< number of fields read from the diag file
 type(esmf_fieldbundle), public        :: input_diag_bundle
                                          !< bundle to hold input diag fields
 type(esmf_fieldbundle), public        :: target_diag_bundle
                                          !< bundle to hold target diag fields

 integer, public                       :: n_hist_fields_2d_patch
                                          !< number of 2d fields read from the hist file
                                          !< to use with patch regridding
 integer, public                       :: n_hist_fields_2d_cons
                                          !< number of 2d fields read from the hist file
                                          !< to use with conservative regridding
 integer, public                       :: n_hist_fields_2d_nstd
                                          !< number of 2d fields read from the hist file
                                          !< to use with nearest source to destination
                                          !< regridding
 integer, public                       :: n_hist_fields_3d_vert
                                          !< number of 3d fields read from the hist file
                                          !< to be defined on input grid vertices and 
                                          !< regridded with bilinear interpolation to 
                                          !< target grid cells 
 integer, public                       :: n_hist_fields_3d_nz
                                          !< number of 3d fields read from the hist file
                                          !< with vertical dimension nVertLevels
 integer, public                       :: n_hist_fields_3d_nzp1
                                          !< number of 3d fields read from the hist file
                                          !< with vertical dimension nVertLevelsp1
 integer, public                       :: n_hist_fields_soil
                                          !< number of soil fields read from the hist file
 real, public                          :: diag_out_interval
                                          !< output_interval from diag file
 integer, public                       :: do_u_interp           
                                          !< whether 3d u is requested for interpolation
 integer, public                       :: do_v_interp           
                                          !< whether 3d v is requested for interpolation
 integer, public                       :: do_u10_interp
                                          !< whether 10-m u is requested for interpolation
 integer, public                       :: do_v10_interp
                                          !< whether 10-m v is requested for interpolation
 integer, public                       :: u10_ind
                                          !< index of u10 in target_diag_bundle
 integer, public                       :: v10_ind
                                          !< index of v10 in target_diag_bundle
 character(50), allocatable, public    :: target_diag_names(:), &
                                          target_hist_names_2d_cons(:), &
                                          target_hist_names_2d_nstd(:), &
                                          target_hist_names_2d_patch(:), &
                                          target_hist_names_3d_nzp1(:), &
                                          target_hist_names_3d_nz(:), &
                                          target_hist_names_3d_vert(:), &
                                          target_hist_names_soil(:)
                                          !< Arrays to hold target field names
 character(50), allocatable, public    :: target_diag_units(:), &
                                          target_hist_units_2d_cons(:), &
                                          target_hist_units_2d_nstd(:), &
                                          target_hist_units_2d_patch(:), &
                                          target_hist_units_3d_nzp1(:), &
                                          target_hist_units_3d_nz(:), &
                                          target_hist_units_3d_vert(:), &
                                          target_hist_units_soil(:)
                                          !< Arrays to hold target field units
 character(200), allocatable, public   :: target_diag_longname(:), &
                                          target_hist_longname_2d_cons(:), &
                                          target_hist_longname_2d_nstd(:), &
                                          target_hist_longname_2d_patch(:), &
                                          target_hist_longname_3d_nzp1(:), &
                                          target_hist_longname_3d_nz(:), &
                                          target_hist_longname_3d_vert(:), &
                                          target_hist_longname_soil(:)
                                          !< Arrays to hold target field longname
 type(esmf_fieldbundle), public        :: input_hist_bundle_2d_patch, &
                                          input_hist_bundle_2d_cons, &
                                          input_hist_bundle_2d_nstd, &
                                          input_hist_bundle_3d_nz, &
                                          input_hist_bundle_3d_nzp1, &
                                          input_hist_bundle_3d_vert, &
                                          input_hist_bundle_soil
                                          !< bundles to hold input hist fields
 type(esmf_fieldbundle), public        :: target_hist_bundle_2d_patch, &
                                          target_hist_bundle_2d_cons, &
                                          target_hist_bundle_2d_nstd, &
                                          target_hist_bundle_3d_nz, &
                                          target_hist_bundle_3d_nzp1, &
                                          target_hist_bundle_3d_vert, &
                                          target_hist_bundle_soil
                                          !< bundles to hold target hist fields

 public :: define_target_grid
 public :: define_input_grid
 public :: cleanup_input_target_grid_data
 contains
 

!> Define input grid object for MPAS input data.
!!
!! @param [in] localpet ESMF local persistent execution thread
!! @param [in] npets  Number of persistent execution threads
!! @author  Larissa Reames CIWRO/NOAA/NSSL/FRDD

 subroutine define_input_grid(localpet,npets)

 use mpi
 use netcdf
 use program_setup, only       : grid_file_input_grid
 implicit none

 character(len=500)           :: the_file, dimname

 integer, intent(in)          :: localpet, npets

 integer                      :: error, i, j, k, rc, n, lmi(1), lma(1)
 integer                      :: NX, NY, format

 integer                               :: ncid,id_var, id_dim, dimsize, nVertThis
 integer                               :: nCells, nVertices, maxEdges, dimids(2)
 integer                               :: cell_start, cell_end, dims(3), temp1, &
                                          temp2, temp3, temp(1), clb(1), cub(1)
 integer, allocatable                  :: elemTypes2(:), vertOnCell(:,:), &
                                          nodesPET(:), nodeIDs_temp(:), &
                                          elementConn_temp(:), elementConn(:), &
                                          unique_nodes(:),node_pets(:),myCells(:)
 real(esmf_kind_r8), allocatable       :: latCell(:), lonCell(:), &
                                          latVert(:), lonVert(:), &
                                          nodeCoords(:), &
                                          nodeCoords_temp(:), &
                                          elemCoords(:), hgt(:)
 real(esmf_kind_r8), pointer           :: data_1d(:), data_1d2(:)
 real(esmf_kind_r8), parameter         :: PI=4.D0*DATAN(1.D0)
 integer :: myCells_num

 the_file = grid_file_input_grid

 if (localpet==0) print*,'- OPEN MPAS INPUT FILE: ',trim(the_file)
!error=nf90_open_par(trim(the_file),NF90_NOWRITE,MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
 error=nf90_open(trim(the_file),nf90_nowrite,ncid) ! CSS
 if (error /=0) call error_handler("OPENING MPAS INPUT FILE",error)

 !Get nCells size
 if (localpet==0) print*,'- READ nCells'
 error = nf90_inq_dimid(ncid,'nCells', id_dim)
 call netcdf_err(error, 'reading nCells id')

 error=nf90_inquire_dimension(ncid,id_dim,len=nCells)
 call netcdf_err(error, 'reading nCells')

 nCells_input = nCells

 !Get nVertices size
  if (localpet==0) print*,'- READ nVertices'
 error = nf90_inq_dimid(ncid,'nVertices',id_dim)
 call netcdf_err(error, 'reading nVertices id')

 error=nf90_inquire_dimension(ncid,id_dim,len=nVertices)
 call netcdf_err(error, 'reading nVertices')

 nVert_input = nVertices
 !Get nVertLevels size
 if (localpet==0) print*,'- READ nVertLevels'
 error = nf90_inq_dimid(ncid,'nVertLevels',id_dim)
 call netcdf_err(error, 'reading nVertLevels id')

 error=nf90_inquire_dimension(ncid,id_dim,len=nz_input)
 call netcdf_err(error, 'reading nVertLevels')

 !Get nVertLevelsP1 size
  if (localpet==0) print*,'- READ nVertLevelsP1'
 error = nf90_inq_dimid(ncid,'nVertLevelsP1',id_dim)
 call netcdf_err(error, 'reading nVertLevelsP1 id')

 error=nf90_inquire_dimension(ncid,id_dim,len=nzp1_input)
 call netcdf_err(error, 'reading nVertLevelsP1')

 !Get maxEdges size
  if (localpet==0) print*,'- READ maxEdges'
 error = nf90_inq_dimid(ncid,'maxEdges',id_dim)
 call netcdf_err(error, 'reading maxEdges id')

 error=nf90_inquire_dimension(ncid,id_dim,len=maxEdges)
 call netcdf_err(error, 'reading maxEdges')

 !Get nSoilLevels size
  if (localpet==0) print*,'- READ nSoilLevels'
 error = nf90_inq_dimid(ncid,'nSoilLevels',id_dim)
 call netcdf_err(error, 'reading nSoilLevels id')

 error=nf90_inquire_dimension(ncid,id_dim,len=nsoil_input)
 call netcdf_err(error, 'reading nSoilLevels')
 
 allocate(latCell(nCells))
 allocate(lonCell(nCells))
 allocate(hgt(nCells))

 allocate(latVert(nVertices))
 allocate(lonVert(nVertices))

 
 allocate(vertOnCell(maxEdges,nCells))
 allocate(zs_target_grid(nsoil_input,1))

 ! GET CELL CENTER LAT/LON
 if (localpet==0) print*,'- READ LONCELL ID'
 error=nf90_inq_varid(ncid, 'lonCell', id_var)
 call netcdf_err(error, 'reading lonCell id')

 if (localpet==0) print*,'- READ LONCELL'
 error=nf90_get_var(ncid, id_var, start=(/1/), count=(/nCells/),values=lonCell)
 call netcdf_err(error, 'reading lonCell')

 if (localpet==0) print*,'- READ LATCELL ID'
 error=nf90_inq_varid(ncid, 'latCell', id_var)
 call netcdf_err(error, 'reading latCell id')

 if (localpet==0) print*,'- READ LATCELL'
  error=nf90_get_var(ncid, id_var, start=(/1/), count=(/nCells/),values=latCell)
 call netcdf_err(error, 'reading latCell')

 ! GET VERTEX LAT/LON
 if (localpet==0) print*,'- READ LONVERTEX ID'
 error=nf90_inq_varid(ncid, 'lonVertex', id_var)
 call netcdf_err(error, 'reading lonVertex id')

 if (localpet==0) print*,'- READ LONVERTEX'
 error=nf90_get_var(ncid, id_var, start=(/1/),count=(/nVertices/),values=lonVert)
 call netcdf_err(error, 'reading lonVertex')
 
 if (localpet==0) print*,'- READ LATVERTEX ID'
 error=nf90_inq_varid(ncid, 'latVertex', id_var)
 call netcdf_err(error, 'reading latVertex id')

 if (localpet==0) print*,'- READ LATVERTEX'
 error=nf90_get_var(ncid, id_var,  start=(/1/),count=(/nVertices/),values=latVert)
 call netcdf_err(error, 'reading latVertex')

  ! SOIL CENTER DEPTHS
 if (localpet==0) print*,'- READ ZS ID'
 error=nf90_inq_varid(ncid, 'zs', id_var)
 call netcdf_err(error, 'reading zs id')

 if (localpet==0) print*,'- READ ZS'
 error=nf90_get_var(ncid, id_var, start=(/1,1/),count=(/nsoil_input,1/),values=zs_target_grid)
 call netcdf_err(error, 'reading ZS')

 if (localpet==0) print*,'- READ HGT'
 error=nf90_inq_varid(ncid, 'ter', id_var)
 call netcdf_err(error, 'reading ter id')
 error=nf90_get_var(ncid, id_var, start=(/1/),count=(/nCells/), values=hgt)
 call netcdf_err(error, 'reading ter')

 if (localpet==0) print*,"- NUMBER OF CELLS ON INPUT GRID ", nCells_input
 if (localpet==0) print*,"- NUMBER OF NODES ON INPUT GRID ", nVert_input

!-----------------------------------------------------------------------
! Create ESMF grid object for the model grid.
!-----------------------------------------------------------------------

 if (localpet==0) print*,'- READ verticesOnCell ID'
 error=nf90_inq_varid(ncid, 'verticesOnCell', id_var)
 call netcdf_err(error, 'reading verticesOnCell id')

 error = nf90_inquire_variable(ncid,id_var,dimids=dims)
 call netcdf_err(error, 'reading verticesOnCell dims')

 if (localpet==0) print*,'- READ verticesOnCell'
 error=nf90_get_var(ncid, id_var, start=(/1,1/),count=(/maxEdges,nCells/),values=vertOnCell)
 call netcdf_err(error, 'reading verticesOnCell')

 error = nf90_close(ncid)

 nVertThis = 0
 k = 0
 if (block_decomp_file=='NULL') then
     !nCellsPerPET = ceiling(real(nCells)/real(npets))
     !cell_start = localpet*nCellsPerPET+1
     !cell_end = min(localpet*nCellsPerPET+nCellsPerPET,nCells)
     call para_range(1, nCells, npets, localpet, cell_start, cell_end)
     nCellsPerPET = cell_end - cell_start + 1
     allocate(elemIDs(nCellsPerPET))
!$OMP PARALLEL DO $PRIVATE(i,j)
     do i = cell_start, cell_end
         j = i - cell_start + 1
         elemIDs(j) = i
     enddo
!$OMP END PARALLEL DO
 else
     call read_block_decomp_file(localpet,npets,block_decomp_file,nCells,elemIDs,nCellsPerPET)
 endif

 allocate(elementConn_temp(maxEdges*nCellsPerPET))
 allocate(elemCoords(2*nCellsPerPET))
 allocate(elemTypes2(nCellsPerPET))
 
 if (localpet==0) print*, "- Create PET-local element connectivity "
!$OMP PARALLEL DO $PRIVATE(i,j,n)
 do i = 1,nCellsPerPET
    j = elemIDs(i)
    elemTypes2(i) = count(vertOnCell(:,j)/=0)
    nVertThis = nVertThis + elemTypes2(i)
    elemCoords(2*i-1) = lonCell(j)*180.0_esmf_kind_r8/PI
    if (elemCoords(2*i-1) > 180.0_esmf_kind_r8) then
        elemCoords(2*i-1) = elemCoords(2*i-1) - 360.0_esmf_kind_r8
    endif
    elemCoords(i*2) = latCell(j)*180.0_esmf_kind_r8/PI
    elementConn_temp(maxEdges*(i-1)+1:maxEdges*i) = vertOnCell(:,j)
 enddo
!$OMP END PARALLEL DO
 call unique_sort(elementConn_temp,maxEdges*nCellsPerPET,nodeIDs)
 nNodesPerPET=size(nodeIDs)
 allocate(nodeCoords(2*nNodesPerPET))
!$OMP PARALLEL DO 
do j = 1,nNodesPerPET
   i = nodeIDs(j)
   nodeCoords(2*j-1) = lonVert(i)*180.0_esmf_kind_r8/PI
   if (nodeCoords(2*j-1) > 180.0_esmf_kind_r8) then
     nodeCoords(2*j-1) = nodeCoords(2*j-1) - 360.0_esmf_kind_r8
   endif
   nodeCoords(2*j) = latVert(i)*180.0_esmf_kind_r8/PI
enddo
!$OMP END PARALLEL DO
allocate(elementConn(nVertThis))
nVertThis = 0
!$OMP PARALLEL DO
do i = 1,nCellsPerPET
  k = elemIDs(i)
  do j = 1,maxEdges
!    if(vertOnCell(j,i)/=0) then
!    temp = FINDLOC(nodeIDs,vertOnCell(j,elemIDs(i)))
   if(vertOnCell(j,k)/=0) then
     temp = FINDLOC(nodeIDs,vertOnCell(j,k))
     elementConn(nVertThis+1) = temp(1)
     nVertThis = nVertThis + 1
   endif
 enddo
enddo
!$OMP END PARALLEL DO
  if (localpet==0) print*, "- CREATE MESH -"
 input_grid = ESMF_MeshCreate(parametricDim=2, &
                     spatialDim=2, &
                     nodeIDs= nodeIDs, &
                     nodeCoords = nodeCoords, &
                     elementIDs = elemIDs, &
                     elementTypes=elemTypes2, &
                     elementConn = elementConn, &
                     elementCoords=elemCoords, &
                     coordSys=ESMF_COORDSYS_SPH_DEG, &
                     rc=rc)
 if (ESMF_LogFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN MeshCreate", rc)


 ! After the creation we are through with the arrays, so they may be deallocated.
 deallocate(elemCoords,elemTypes2)
 deallocate(nodeCoords)
 deallocate(elementConn_temp, elementConn)
 !deallocate(nodeIDs_temp)


 !-----------------------------------------------------------------------
 ! Create lat/lon arrays on input element grid
 !-----------------------------------------------------------------------

 if (localpet==0) print*,"- CALL FieldCreate FOR INPUT GRID CELL LATITUDE."
 cell_latitude_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   meshloc=ESMF_MESHLOC_ELEMENT, &
                                   name="input_grid_cell_latitude", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", error)

 if (localpet==0) print*,"- CALL FieldCreate FOR INPUT GRID CELL LONGITUDE."
 cell_longitude_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   meshloc=ESMF_MESHLOC_ELEMENT, &
                                   name="input_grid_cell_longitude", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", error)

 call ESMF_FieldGet(cell_latitude_input_grid, farrayPtr=data_1d,rc = rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", error)

 call ESMF_FieldGet(cell_longitude_input_grid,farrayPtr=data_1d2, &
                    computationalLBound=clb, computationalUBound=cub, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldGet", error)

 
 do j = 1, nCellsPerPET
    i = elemIDs(j)
    data_1d(j) = latCell(i)
    data_1d2(j) = lonCell(i)
 enddo

 nullify(data_1d,data_1d2)

 !-----------------------------------------------------------------------
 ! Create lat/lon arrays on input node grid
 !-----------------------------------------------------------------------

 if (localpet==0) print*,"- CALL FieldCreate FOR INPUT GRID VERTEX LATITUDE."
 node_latitude_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   meshloc=ESMF_MESHLOC_NODE, &
                                   name="input_grid_node_latitude", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", error)

 if (localpet==0) print*,"- CALL FieldCreate FOR INPUT GRID VERTEX LONGITUDE."
 node_longitude_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   meshloc=ESMF_MESHLOC_NODE, &
                                   name="input_grid_node_longitude", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", error)

 call ESMF_FieldGet(node_latitude_input_grid, farrayPtr=data_1d,rc = rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldGet", error)

 call ESMF_FieldGet(node_longitude_input_grid,farrayPtr=data_1d2, &
                    computationalLBound=clb, computationalUBound=cub, rc=rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldGet", error)

 !call unique_sort(nodeIDs,k,unique_nodes)
 call ESMF_MeshGet(input_grid,numOwnedNodes=n)
 allocate(unique_nodes(n))
 allocate(node_pets(nNodesPerPET))
 call ESMF_MeshGet(input_grid,nodeOwners=node_pets)
 j = 1
 do i = 1,nNodesPerPET
   if (node_pets(i)==localpet) then
     unique_nodes(j) = nodeIDs(i)
     j = j + 1
   endif
 enddo
 
 do j = 1,n
    i = unique_nodes(j)
    data_1d(j) = latVert(i)
    data_1d2(j) = lonVert(i)
 enddo
 nNodesPerPET = n
 nullify(data_1d,data_1d2)

 if (localpet==0) print*,"- CALL FieldCreate FOR INPUT GRID HGT."
 hgt_input_grid = ESMF_FieldCreate(input_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   meshloc=ESMF_MESHLOC_ELEMENT, &
                                   name="input_grid_hgt", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", error)

 call ESMF_FieldGet(hgt_input_grid, farrayPtr=data_1d,rc = rc)
  if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldGet", error)

 
 do j = 1,nCellsPerPET
    i = elemIDs(j)
    data_1d(j) = hgt(i)
 enddo
 nullify(data_1d)
 deallocate(lonCell, latCell, lonVert, latVert,vertOnCell)


 end subroutine define_input_grid

!> Setup the esmf grid object for the target grid.
!!
!! @param [in] localpet ESMF local persistent execution thread
!! @param [in] npets Number of persistent execution threads
!! @author Larissa Reames CIWRO/NOAA/NSSL/FRDD
 subroutine define_target_grid(localpet, npets)
 
 use program_setup, only : target_grid_type
 implicit none

 integer, intent(in) :: localpet, npets
 if(trim(target_grid_type) == 'file') then
    call define_target_grid_file(localpet,npets)
 else
    call define_target_grid_params(localpet,npets)
 endif
 
 end subroutine define_target_grid
 
 subroutine define_target_grid_params(localpet,npets)
 
 use netcdf

 use llxy_module

 implicit none

 integer, intent(in)          :: localpet, npets

 integer                      :: error, extra, i, j, clb(2), cub(2)

 real(esmf_kind_r8), allocatable       :: latitude_one(:,:), longitude_one(:,:), &
                                          latitude_corner_one(:,:), &
                                          longitude_corner_one(:,:), &
                                          longitude_u_one(:,:), latitude_u_one(:,:), &
                                          longitude_v_one(:,:), latitude_v_one(:,:), &
                                          mapfac_m_one(:,:), mapfac_u_one(:,:), &
                                          mapfac_v_one(:,:)
                                          
 integer                               :: ncid,id_var, id_dim
 real(esmf_kind_r8), pointer           :: lat_src_ptr(:,:), lon_src_ptr(:,:),& 
                                           mapptr(:,:), cosa_ptr(:,:), sina_ptr(:,:)

 type(esmf_polekind_flag)              :: polekindflag(2)
 
 ip1_target = i_target + 1
 jp1_target = j_target + 1
 
 !----------------------------------------------------------------------
 ! Fill proj object for target grid projection
 !----------------------------------------------------------------------
 call push_source_projection(proj_code, stand_lon, truelat1, truelat2, &
              dxkm, dykm, dlatdeg, dlondeg, known_x, known_y, &
              known_lat, known_lon)
                              
!-----------------------------------------------------------------------
! Create ESMF grid object for the model grid.
!-----------------------------------------------------------------------

 if (.not. is_regional ) then
     polekindflag(1:2) = ESMF_POLEKIND_MONOPOLE
     print*,"- CALL GridCreate1PeriDim FOR INPUT GRID."
     target_grid = ESMF_GridCreate1PeriDim(minIndex=(/1,1/), &
            maxIndex=(/i_target,j_target/), &
            polekindflag=polekindflag, &
            periodicDim=1, &
            poleDim=2,  &
            coordSys=ESMF_COORDSYS_SPH_DEG, &
            regDecomp=(/1,npets/),  &
            indexflag=ESMF_INDEX_GLOBAL, rc=error)
     if(ESMF_logFoundError(rcToCheck=error, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
       call error_handler("IN GridCreate1PeriDim", error)
 else
     if (localpet==0) print*,"- CALL GridCreateNoPeriDim FOR TARGET MODEL GRID"
     target_grid = ESMF_GridCreateNoPeriDim(maxIndex=(/i_target,j_target/), &
               indexflag=ESMF_INDEX_GLOBAL, &
               rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN GridCreateNoPeriDim", error)
 endif
    
 if (localpet==0) print*,"- CALL GridAddCoord FOR INPUT GRID CENTER."
 call ESMF_GridAddCoord(target_grid, &
                        staggerloc=ESMF_STAGGERLOC_CENTER, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridAddCoord", error)
    
  if (localpet==0) print*,"- CALL GridAddCoord FOR INPUT GRID CORNER."
 call ESMF_GridAddCoord(target_grid, &
                        staggerloc=ESMF_STAGGERLOC_CORNER, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridAddCoord", error)
    
  if (localpet==0) print*,"- CALL GridAddCoord FOR INPUT GRID EDGE1."
 call ESMF_GridAddCoord(target_grid, &
                        staggerloc=ESMF_STAGGERLOC_EDGE1, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridAddCoord", error)
    
  if (localpet==0) print*,"- CALL GridAddCoord FOR INPUT GRID EDGE2."
 call ESMF_GridAddCoord(target_grid, &
                        staggerloc=ESMF_STAGGERLOC_EDGE2, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridAddCoord", error)
                              
   
!-----------------------------------------------------------------------
! Generate lat/lon array values for various staggers
!----------------------------------------------------------------------- 

 !Grid centers
 call ESMF_GridGetCoord(target_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CENTER, &
                        coordDim=2, &
                        computationalLBound=clb, &
                        computationalUBound=cub, &
                        farrayPtr=lat_src_ptr, &
                        rc=error)  
                          
 allocate(latitude_one(clb(1):cub(1),clb(2):cub(2)))
 allocate(longitude_one(clb(1):cub(1),clb(2):cub(2)))
 allocate(mapfac_m_one(clb(1):cub(1),clb(2):cub(2)))                   
 call get_lat_lon_fields(latitude_one, longitude_one, clb(1),clb(2),cub(1),cub(2),M) 
 call get_map_factor(latitude_one, longitude_one, mapfac_m_one, mapfac_m_one, clb(1),clb(2),cub(1),cub(2))
 
 !y-dir stagger (V)
 
  call ESMF_GridGetCoord(target_grid, &
                          staggerLoc=ESMF_STAGGERLOC_EDGE2, &
                          coordDim=2, &
                          computationalLBound=clb, &
                          computationalUBound=cub, &
                          farrayPtr=lat_src_ptr, &
                          rc=error)  
                          
 allocate(latitude_v_one(clb(1):cub(1),clb(2):cub(2)))
 allocate(longitude_v_one(clb(1):cub(1),clb(2):cub(2)))  
 allocate(mapfac_v_one(clb(1):cub(1),clb(2):cub(2))) 
 call get_lat_lon_fields(latitude_v_one, longitude_v_one, clb(1),clb(2),cub(1),cub(2), V)
 call get_map_factor(latitude_v_one, longitude_v_one, mapfac_v_one, mapfac_v_one, clb(1),clb(2),cub(1),cub(2))
 
 
 !x-dir stagger (U)
 
  call ESMF_GridGetCoord(target_grid, &
                          staggerLoc=ESMF_STAGGERLOC_EDGE1, &
                          coordDim=2, &
                          computationalLBound=clb, &
                          computationalUBound=cub, &
                          farrayPtr=lat_src_ptr, &
                          rc=error)  
                          
 allocate(latitude_u_one(clb(1):cub(1),clb(2):cub(2)))
 allocate(longitude_u_one(clb(1):cub(1),clb(2):cub(2)))
 allocate(mapfac_u_one(clb(1):cub(1),clb(2):cub(2))) 
 call get_lat_lon_fields(latitude_u_one, longitude_u_one, clb(1),clb(2),cub(1),cub(2), U)
 call get_map_factor(latitude_u_one, longitude_u_one, mapfac_u_one, mapfac_u_one, clb(1),clb(2),cub(1),cub(2))
  
 !corner
call ESMF_GridGetCoord(target_grid, &
                          staggerLoc=ESMF_STAGGERLOC_CORNER, &
                          coordDim=2, &
                          computationalLBound=clb, &
                          computationalUBound=cub, &
                          farrayPtr=lat_src_ptr, &
                          rc=error)
                          
 allocate(latitude_corner_one(clb(1):cub(1),clb(2):cub(2)))
 allocate(longitude_corner_one(clb(1):cub(1),clb(2):cub(2)))
 call get_lat_lon_fields(latitude_corner_one, longitude_corner_one, clb(1),clb(2),cub(1),cub(2),CORNER)
 !print*, minval(latitude_corner_one), maxval(latitude_corner_one), minval(longitude_corner_one), maxval(longitude_corner_one)
!-----------------------------------------------------------------------
! Read the mask and lat/lons.
!-----------------------------------------------------------------------
! Grid Centers
 if (localpet==0) print*,"- CALL FieldCreate FOR TARGET GRID LATITUDE."
 latitude_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="target_grid_latitude", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", error)


 if (localpet==0) print*,"- CALL FieldCreate FOR TARGET GRID LONGITUDE."
 longitude_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="target_grid_longitude", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", error)

 if (localpet==0) print*, "- CALL FieldGet FOR TARGET GRID LATITUDE."
 call ESMF_FieldGet(latitude_target_grid, farrayPtr=lat_src_ptr,rc=error)
  if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldGET", error)

 if (localpet==0) print*, "- CALL FieldGet FOR TARGET GRID LONGITUDE."
 call ESMF_FieldGet(longitude_target_grid, farrayPtr=lon_src_ptr, &
                        computationalLBound=clb, computationalUBound=cub, &
                        rc=error)
  if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldGet", error)

  do j = clb(2),cub(2)
  do i = clb(1),cub(1)
        lat_src_ptr(i,j) = latitude_one(i,j)
        lon_src_ptr(i,j) = longitude_one(i,j)
  enddo
  enddo
  nullify(lat_src_ptr, lon_src_ptr)
! E-W Stagger
if (localpet==0) print*,"- CALL FieldCreate FOR TARGET GRID LATITUDE."
 latitude_u_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE1, &
                                   name="target_grid_latitude", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", error)


 if (localpet==0) print*,"- CALL FieldCreate FOR TARGET GRID LONGITUDE."
 longitude_u_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE1, &
                                   name="target_grid_longitude", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", error)

 if (localpet==0) print*, "- CALL FieldGet FOR TARGET GRID LATITUDE."
 call ESMF_FieldGet(latitude_u_target_grid, farrayPtr=lat_src_ptr,rc=error)
  if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldGET", error)

 if (localpet==0) print*, "- CALL FieldGet FOR TARGET GRID LONGITUDE."
 call ESMF_FieldGet(longitude_u_target_grid, farrayPtr=lon_src_ptr, &
                        computationalLBound=clb, computationalUBound=cub, &
                        rc=error)
  if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldGet", error)
  do j = clb(2),cub(2)
  do i = clb(1),cub(1)
        lat_src_ptr(i,j) = latitude_u_one(i,j)
        lon_src_ptr(i,j) = longitude_u_one(i,j)
  enddo
  enddo

  nullify(lat_src_ptr, lon_src_ptr)
 
 ! N-S Stagger
if (localpet==0) print*,"- CALL FieldCreate FOR TARGET GRID LATITUDE."
 latitude_v_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE2, &
                                   name="target_grid_latitude", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", error)


 if (localpet==0) print*,"- CALL FieldCreate FOR TARGET GRID LONGITUDE."
 longitude_v_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE2, &
                                   name="target_grid_longitude", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", error)

 if (localpet==0) print*, "- CALL FieldGet FOR TARGET GRID LATITUDE."
 call ESMF_FieldGet(latitude_v_target_grid, farrayPtr=lat_src_ptr,rc=error)
  if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldGET", error)

 if (localpet==0) print*, "- CALL FieldGet FOR TARGET GRID LONGITUDE."
 call ESMF_FieldGet(longitude_v_target_grid, farrayPtr=lon_src_ptr, &
                        computationalLBound=clb, computationalUBound=cub, &
                        rc=error)
  if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldGet", error)
  !print*, localpet, cub(1), cub(2)
  do j = clb(2),cub(2)
  do i = clb(1),cub(1)
        lat_src_ptr(i,j) = latitude_v_one(i,j)
        lon_src_ptr(i,j) = longitude_v_one(i,j)
  enddo
  enddo
  nullify(lat_src_ptr, lon_src_ptr)

 if (localpet==0) print*,"- CALL FieldCreate FOR TARGET GRID HGT."
 hgt_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="target_grid_hgt", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", error)

 
 !Create and fill coordinates
 if (localpet==0) print*,"- CALL GridGetCoord FOR INPUT GRID X-COORD."
   nullify(lon_src_ptr)
   call ESMF_GridGetCoord(target_grid, &
                          staggerLoc=ESMF_STAGGERLOC_CENTER, &
                          coordDim=1, &
                          farrayPtr=lon_src_ptr, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN GridGetCoord", error)

   if (localpet==0) print*,"- CALL GridGetCoord FOR INPUT GRID Y-COORD."
   nullify(lat_src_ptr)
   call ESMF_GridGetCoord(target_grid, &
                          staggerLoc=ESMF_STAGGERLOC_CENTER, &
                          coordDim=2, &
                          computationalLBound=clb, &
                          computationalUBound=cub, &
                          farrayPtr=lat_src_ptr, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN GridGetCoord", error)

    do j = clb(2),cub(2)
      do i = clb(1), cub(1)
        lon_src_ptr(i,j)=longitude_one(i,j)
        lat_src_ptr(i,j)=latitude_one(i,j)
      enddo
    enddo

  nullify(lon_src_ptr)
  nullify(lat_src_ptr)

 if (localpet==0) print*,"- CALL GridGetCoord FOR TARGET GRID X-COORD."

 call ESMF_GridGetCoord(target_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CORNER, &
                        coordDim=1, &
                        farrayPtr=lon_src_ptr, rc=error)
 if(ESMF_logFoundError(rcToCheck=error, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridGetCoord", error)

 if (localpet==0) print*,"- CALL GridGetCoord FOR TARGET GRID Y-COORD."

 call ESMF_GridGetCoord(target_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CORNER, &
                        coordDim=2, &
                        computationalLBound=clb, &
                        computationalUBound=cub, &
                        farrayPtr=lat_src_ptr, rc=error)
 if(ESMF_logFoundError(rcToCheck=error, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridGetCoord", error)

 do j = clb(2),cub(2)
      do i = clb(1), cub(1)
          lon_src_ptr(i,j)=longitude_corner_one(i,j)
          lat_src_ptr(i,j)=latitude_corner_one(i,j)
      enddo
 enddo

 nullify(lon_src_ptr,lat_src_ptr)
 if (localpet==0) print*,"- CALL GridGetCoord FOR TARGET GRID EDGE1 X-COORD."
 call ESMF_GridGetCoord(target_grid, &
                        staggerLoc=ESMF_STAGGERLOC_EDGE1, &
                        coordDim=1, &
                        farrayPtr=lon_src_ptr, rc=error)
 if(ESMF_logFoundError(rcToCheck=error, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__))&
    call error_handler("IN GridGetCoord", error)

 if (localpet==0) print*,"- CALL GridGetCoord FOR TARGET GRID EDGE1 Y-COORD."

 call ESMF_GridGetCoord(target_grid, &
                        staggerLoc=ESMF_STAGGERLOC_EDGE1, &
                        coordDim=2, &
                        computationalLBound=clb, &
                        computationalUBound=cub, &
                        farrayPtr=lat_src_ptr, rc=error)
 if(ESMF_logFoundError(rcToCheck=error, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__))&
    call error_handler("IN GridGetCoord", error)
 
 do j = clb(2),cub(2)
      do i = clb(1), cub(1)
          lon_src_ptr(i,j)=longitude_u_one(i,j)
          lat_src_ptr(i,j)=latitude_u_one(i,j)
      enddo
 enddo
 
 nullify(lon_src_ptr,lat_src_ptr)
if (localpet==0) print*,"- CALL GridGetCoord FOR TARGET GRID EDGE2 X-COORD."
 call ESMF_GridGetCoord(target_grid, &
                        staggerLoc=ESMF_STAGGERLOC_EDGE2, &
                        coordDim=1, &
                        farrayPtr=lon_src_ptr, rc=error)
 if(ESMF_logFoundError(rcToCheck=error, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__))&
    call error_handler("IN GridGetCoord", error)

 if (localpet==0) print*,"- CALL GridGetCoord FOR TARGET GRID EDGE2 Y-COORD."

 call ESMF_GridGetCoord(target_grid, &
                        staggerLoc=ESMF_STAGGERLOC_EDGE2, &
                        coordDim=2, &
                        computationalLBound=clb, &
                        computationalUBound=cub, &
                        farrayPtr=lat_src_ptr, rc=error)
 if(ESMF_logFoundError(rcToCheck=error, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__))&
    call error_handler("IN GridGetCoord", error)

 do j = clb(2),cub(2)
      do i = clb(1), cub(1)
          lon_src_ptr(i,j)=longitude_v_one(i,j)
          lat_src_ptr(i,j)=latitude_v_one(i,j)
      enddo
 enddo

if (localpet==0) print*,"- CALL FieldCreate FOR TARGET GRID mapfac_m."
 mapfac_m_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="target_grid_hgt", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", error)

 if (localpet==0) print*,"- CALL FieldGet FOR TARGET GRID mapfac_m. "
 call ESMF_FieldGet(mapfac_m_target_grid, farrayPtr=mapptr, &
                    computationalLBound=clb,computationalUBound=cub, &
                    rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
   call error_handler("IN FieldGet", error)

 do j=clb(2),cub(2)
 do i=clb(1),cub(1)
        mapptr(i,j) = mapfac_m_one(i,j)
 enddo
 enddo

 if (localpet==0) print*,"- CALL FieldCreate FOR TARGET GRID mapfac_u."
 mapfac_u_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE1, &
                                   name="target_grid_hgt", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", error)

 if (localpet==0) print*,"- CALL FieldGet FOR TARGET GRID mapfac_u. "
 call ESMF_FieldGet(mapfac_u_target_grid, farrayPtr=mapptr, &
                    computationalLBound=clb,computationalUBound=cub, &
                     rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
   call error_handler("IN FieldGet", error)

 do j=clb(2),cub(2)
 do i=clb(1),cub(1)
        mapptr(i,j) = mapfac_u_one(i,j)
 enddo
 enddo

 if (localpet==0) print*,"- CALL FieldCreate FOR TARGET GRID mapfac_v."
 mapfac_v_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE2, &
                                   name="target_grid_hgt", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", error)

 if (localpet==0) print*,"- CALL FieldGet FOR TARGET GRID mapfac_v. "
 call ESMF_FieldGet(mapfac_v_target_grid, farrayPtr=mapptr, &
                    computationalLBound=clb,computationalUBound=cub, &
                    rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
   call error_handler("IN FieldGet", error)

 do j=clb(2),cub(2)
 do i=clb(1),cub(1)
        mapptr(i,j) = mapfac_v_one(i,j)
 enddo
 enddo
 
 ! Define some extra projection-related values 
 call xytoll(i_target/2.0,j_target/2.0,ref_lat,ref_lon,M)

 ! ------------------------------------------------------------------------
 ! Create wind rotation angle arrays for lambert conformal grids
 ! ------------------------------------------------------------------------

 if (proj_code==PROJ_LC) then
    cosa_target_grid = ESMF_FieldCreate(target_grid, &
                                        typekind=ESMF_TYPEKIND_R8, &
                                        staggerloc=ESMF_STAGGERLOC_CENTER, &
                                        name="cosalpha", &
                                        rc=error)
    sina_target_grid = ESMF_FieldCreate(target_grid, &
                                        typekind=ESMF_TYPEKIND_R8, &
                                        staggerloc=ESMF_STAGGERLOC_CENTER, &
                                        name="sinalpha", &
                                        rc=error)
    cosa_u_target_grid = ESMF_FieldCreate(target_grid, &
                                        typekind=ESMF_TYPEKIND_R8, &
                                        staggerloc=ESMF_STAGGERLOC_EDGE1, &
                                        name="cosalpha_u", &
                                        rc=error)
    sina_u_target_grid = ESMF_FieldCreate(target_grid, &
                                        typekind=ESMF_TYPEKIND_R8, &
                                        staggerloc=ESMF_STAGGERLOC_EDGE1, &
                                        name="sinalpha_u", &
                                        rc=error)
    cosa_v_target_grid = ESMF_FieldCreate(target_grid, &
                                        typekind=ESMF_TYPEKIND_R8, &
                                        staggerloc=ESMF_STAGGERLOC_EDGE2, &
                                        name="cosalpha_v", &
                                        rc=error)
    sina_v_target_grid = ESMF_FieldCreate(target_grid, &
                                        typekind=ESMF_TYPEKIND_R8, &
                                        staggerloc=ESMF_STAGGERLOC_EDGE2, &
                                        name="sinalpha_v", &
                                        rc=error)
    ! Set angle on grid center
    call ESMF_FieldGet(sina_target_grid, farrayPtr=sina_ptr, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldGet", error)
    call ESMF_FieldGet(cosa_target_grid, farrayPtr=cosa_ptr, &
                        computationalLBound=clb,computationalUBound=cub, &
                        rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldGet", error)

    call get_rotang(latitude_one, longitude_one, cosa_ptr, sina_ptr, &
                         clb(1), clb(2), cub(1), cub(2)) 
    nullify(sina_ptr,cosa_ptr)

    ! Set angle on grid edge1
    call ESMF_FieldGet(sina_u_target_grid, farrayPtr=sina_ptr, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldGet", error)
    call ESMF_FieldGet(cosa_u_target_grid, farrayPtr=cosa_ptr, &
                        computationalLBound=clb,computationalUBound=cub, &
                        rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldGet", error)

    call get_rotang(latitude_u_one, longitude_u_one, cosa_ptr, sina_ptr, &
                         clb(1), clb(2), cub(1), cub(2))
    nullify(cosa_ptr,sina_ptr)

    ! Set angle on grid edge2
    call ESMF_FieldGet(sina_v_target_grid, farrayPtr=sina_ptr, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldGet", error)
    call ESMF_FieldGet(cosa_v_target_grid, farrayPtr=cosa_ptr, &
                        computationalLBound=clb,computationalUBound=cub, &
                        rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldGet", error)

    call get_rotang(latitude_v_one, longitude_v_one, cosa_ptr, sina_ptr, &
                         clb(1), clb(2), cub(1), cub(2))
    nullify(cosa_ptr,sina_ptr)
 endif

 nullify(lon_src_ptr)
 nullify(lat_src_ptr)
 deallocate(longitude_one)
 deallocate(latitude_one)
 deallocate(longitude_corner_one)
 deallocate(latitude_corner_one)
 deallocate(latitude_u_one)
 deallocate(longitude_u_one)
 deallocate(latitude_v_one)
 deallocate(longitude_v_one)
 deallocate(mapfac_m_one)
 deallocate(mapfac_u_one)
 deallocate(mapfac_v_one)
 
 end subroutine define_target_grid_params
 
 subroutine define_target_grid_file(localpet,npets)

 use netcdf
 use mpi
 use program_setup, only       : file_target_grid

 implicit none

 character(len=500)           :: the_file

 integer, intent(in)          :: localpet, npets

 integer                      :: error, extra, i, j, clb(2), cub(2), &
                                 starts(2), counts(2)

 real(esmf_kind_r8), allocatable       :: latitude(:,:), longitude(:,:), &
                                          dum2d(:,:), templat(:,:), templon(:,:), &
                                          mapfac_temp(:,:)
 integer                               :: ncid,id_var, id_dim
 real(esmf_kind_r8), pointer           :: lat_src_ptr(:,:), lon_src_ptr(:,:), &
                                          mapptr(:,:), hgtptr(:,:), &
                                          clat_src_ptr(:,:), clon_src_ptr(:,:)


 the_file = file_target_grid

 if (localpet==0) print*,'- OPEN WRF INPUT FILE: ',trim(the_file)
!error=nf90_open_par(trim(the_file),NF90_NOWRITE,MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
 error=nf90_open(trim(the_file),nf90_nowrite,ncid) ! CSS

 if (error /=0) call error_handler("OPENING WRF INPUT FILE",error)

 if (localpet==0) print*,'- READ WEST_EAST ID'
 error=nf90_inq_dimid(ncid, 'west_east', id_dim)
 call netcdf_err(error, 'reading west_east id')

 if (localpet==0) print*,'- READ WEST_EAST'
 error=nf90_inquire_dimension(ncid,id_dim,len=i_target)
 call netcdf_err(error, 'reading west_east')

 if (localpet==0) print*,'- READ SOUTH_NORTH ID'
 error=nf90_inq_dimid(ncid, 'south_north', id_dim)
 call netcdf_err(error, 'reading south_north id')

 if (localpet==0) print*,'- READ SOUTH_NORTH'
 error=nf90_inquire_dimension(ncid,id_dim,len=j_target)
 call netcdf_err(error, 'reading south_north')

 if (localpet==0) print*,'- READ GLOBAL ATTRIBUTE DX'
 error = nf90_get_att(ncid,NF90_GLOBAL,'DX',dx)
 call netcdf_err(error, 'reading dx')

 if (localpet==0) print*,"- I/J DIMENSIONS OF THE TARGET GRID TILES ", i_target, j_target

 ip1_target = i_target + 1
 jp1_target = j_target + 1

 if (localpet==0) print*, '- READING GLOBAL ATTRIBUTES'
 error = nf90_get_att(ncid, NF90_GLOBAL, 'CEN_LAT', ref_lat)
 call netcdf_err(error, 'GETTING CEN_LAT GLOBAL ATTRIBUTE')

 error = nf90_get_att(ncid, NF90_GLOBAL, 'CEN_LON', ref_lon)
 call netcdf_err(error, 'GETTING CEN_LON GLOBAL ATTRIBUTE')

 error = nf90_get_att(ncid, NF90_GLOBAL, 'TRUELAT1', truelat1)
 call netcdf_err(error, 'GETTING TRUELAT1 GLOBAL ATTRIBUTE')

 error = nf90_get_att(ncid, NF90_GLOBAL, 'TRUELAT2', truelat2)
 call netcdf_err(error, 'GETTING TRUELAT2 GLOBAL ATTRIBUTE')

 error = nf90_get_att(ncid, NF90_GLOBAL, 'MOAD_CEN_LAT', ref_lat)
 call netcdf_err(error, 'GETTING MOAD_CEN_LAT GLOBAL ATTRIBUTE')

 error = nf90_get_att(ncid, NF90_GLOBAL, 'STAND_LON', stand_lon)
 call netcdf_err(error, 'GETTING STAND_LON GLOBAL ATTRIBUTE')

 error = nf90_get_att(ncid, NF90_GLOBAL, 'POLE_LAT', pole_lat)
 call netcdf_err(error, 'GETTING POLELAT GLOBAL ATTRIBUTE')

 error = nf90_get_att(ncid, NF90_GLOBAL, 'POLE_LON', pole_lon)
 call netcdf_err(error, 'GETTING POLE_LON GLOBAL ATTRIBUTE')

 error = nf90_get_att(ncid, NF90_GLOBAL, 'MAP_PROJ', proj_code)
 call netcdf_err(error, 'GETTING MAP_PROJ GLOBAL ATTRIBUTE')

 error = nf90_get_att(ncid, NF90_GLOBAL, 'MAP_PROJ_CHAR', map_proj_char)
 call netcdf_err(error, 'GETTING MAP_PROJ GLOBAL ATTRIBUTE')
 if (error .ne. 0) then
   if (proj_code == 1) then
     map_proj_char = "Lambert Conformal"
   else
     map_proj_char = "Lat/Lon"
   endif
 endif
 
!-----------------------------------------------------------------------
! Create ESMF grid object for the model grid.
!-----------------------------------------------------------------------

 if (localpet==0) print*,"- CALL GridCreateNoPeriDim FOR TARGET MODEL GRID"
 target_grid = ESMF_GridCreateNoPeriDim(maxIndex=(/i_target,j_target/), &
                                       indexflag=ESMF_INDEX_GLOBAL, &
                                       rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridCreateNoPeriDim", error)


!-----------------------------------------------------------------------
! Create needed grid coordinates
!-----------------------------------------------------------------------
    
  if (localpet==0) print*,"- CALL GridAddCoord FOR INPUT GRID CENTER."
 call ESMF_GridAddCoord(target_grid, &
                        staggerloc=ESMF_STAGGERLOC_CENTER, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridAddCoord", error)
    
  if (localpet==0) print*,"- CALL GridAddCoord FOR INPUT GRID CORNER."
 call ESMF_GridAddCoord(target_grid, &
                        staggerloc=ESMF_STAGGERLOC_CORNER, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridAddCoord", error)
    
  if (localpet==0) print*,"- CALL GridAddCoord FOR INPUT GRID EDGE1."
 call ESMF_GridAddCoord(target_grid, &
                        staggerloc=ESMF_STAGGERLOC_EDGE1, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridAddCoord", error)
    
  if (localpet==0) print*,"- CALL GridAddCoord FOR INPUT GRID EDGE2."
 call ESMF_GridAddCoord(target_grid, &
                        staggerloc=ESMF_STAGGERLOC_EDGE2, rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridAddCoord", error)
    
!----------------------------------------------------------
!  Read in coordinate values and set on grid
!----------------------------------------------------------
    
!------- Grid center coordinates
    
  if (localpet==0) print*,"- CALL GridGetCoord FOR INPUT GRID X-COORD."
  nullify(lon_src_ptr)
  call ESMF_GridGetCoord(target_grid, &
                      staggerLoc=ESMF_STAGGERLOC_CENTER, &
                      coordDim=1, &
                      farrayPtr=lon_src_ptr, rc=error)
  if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
   call error_handler("IN GridGetCoord", error)

  if (localpet==0) print*,"- CALL GridGetCoord FOR INPUT GRID Y-COORD."
  nullify(lat_src_ptr)
  call ESMF_GridGetCoord(target_grid, &
                      staggerLoc=ESMF_STAGGERLOC_CENTER, &
                      coordDim=2, &
                      computationalLBound=clb, &
                      computationalUBound=cub, &
                      farrayPtr=lat_src_ptr, rc=error)
if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
  call error_handler("IN GridGetCoord", error)

 allocate(templat(clb(1):cub(1),clb(2):cub(2)))
 allocate(templon(clb(1):cub(1),clb(2):cub(2)))
 starts = (/clb(1),clb(2)/)
 counts = (/cub(1)-clb(1)+1,cub(2)-clb(2)+1/)
 
 if (localpet==0) print*,'- READ LONGITUDE ID'
 error=nf90_inq_varid(ncid, 'XLONG', id_var)
 if (error /= NF90_NOERR) then
   error=nf90_inq_varid(ncid, 'XLONG_M', id_var)
   call netcdf_err(error, 'reading longitude id')
 endif

 if (localpet==0) print*,'- READ LONGITUDE'
 error=nf90_get_var(ncid, id_var, start=starts,count=counts,values=templon)
 call netcdf_err(error, 'reading longitude')
 
 if (localpet==0) print*,'- READ LATITUDE ID'
 error=nf90_inq_varid(ncid, 'XLAT', id_var)
 if (error .ne. NF90_NOERR) then
   error=nf90_inq_varid(ncid, 'XLAT_M', id_var)
   call netcdf_err(error, 'reading latitude id')
 endif

 if (localpet==0) print*,'- READ LATITUDE'
 error=nf90_get_var(ncid, id_var, start=starts,count=counts,values=templat)
 call netcdf_err(error, 'reading latitude')
    
    do j = clb(2),cub(2)
      do i = clb(1), cub(1)
        lon_src_ptr(i,j)=real(templon(i,j),esmf_kind_r8)
        lat_src_ptr(i,j)=real(templat(i,j),esmf_kind_r8)
      enddo
    enddo

  nullify(lon_src_ptr)
  nullify(lat_src_ptr)
  
  if (localpet==0) print*,"- CALL FieldCreate FOR TARGET GRID LONGITUDE."
 longitude_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="target_grid_longitude", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldCreate", error)

 if (localpet==0) print*,"- CALL FieldCreate FOR TARGET GRID LATITUDE."
 latitude_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="target_grid_latitude", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", error)
    
 call ESMF_FieldGet(latitude_target_grid, farrayptr=lat_src_ptr,rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldGet", error)
 call ESMF_FieldGet(longitude_target_grid, farrayptr=lon_src_ptr,rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldGet", error)
 
 do j = clb(2),cub(2)
 do i = clb(1), cub(1)
    lon_src_ptr(i,j)=real(templon(i,j),esmf_kind_r8)
    lat_src_ptr(i,j)=real(templat(i,j),esmf_kind_r8)
 enddo
 enddo  
  
!---------- Grid corners coordinate creation

  if (localpet==0) print*,"- CALL GridGetCoord FOR INPUT CORNERS GRID X-COORD."
   nullify(lon_src_ptr)
   call ESMF_GridGetCoord(target_grid, &
                          staggerLoc=ESMF_STAGGERLOC_CORNER, &
                          coordDim=1, &
                          farrayPtr=clon_src_ptr, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN GridGetCoord", error)

   if (localpet==0) print*,"- CALL GridGetCoord FOR INPUT CORNERS GRID Y-COORD."
   nullify(lat_src_ptr)
   call ESMF_GridGetCoord(target_grid, &
                          staggerLoc=ESMF_STAGGERLOC_CORNER, &
                          coordDim=2, &
                          computationalLBound=clb, &
                          computationalUBound=cub, &
                          farrayPtr=clat_src_ptr, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN GridGetCoord", error)
      
  allocate(latitude(i_target,j_target))
  allocate(longitude(i_target,j_target))
  call ESMF_FieldGet(latitude_target_grid, farrayPtr=lat_src_ptr, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGet", error)
      

  call ESMF_FieldGet(longitude_target_grid, farrayPtr=lon_src_ptr, rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldGet", error)
   
  call get_cell_corners(lat_src_ptr, lon_src_ptr, clat_src_ptr, clon_src_ptr, dx, clb, cub)
    
  nullify(lon_src_ptr)
  nullify(lat_src_ptr)
  nullify(clon_src_ptr)
  nullify(clat_src_ptr)
  deallocate(templat,templon)
  deallocate(latitude,longitude)
  
! ------- E-W Stagger grid 

  ! Get coordinate pointers from grid
  if (localpet==0) print*,"- CALL GridGetCoord FOR INPUT U GRID X-COORD."
   nullify(lon_src_ptr)
   call ESMF_GridGetCoord(target_grid, &
                          staggerLoc=ESMF_STAGGERLOC_EDGE1, &
                          coordDim=1, &
                          farrayPtr=lon_src_ptr, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN GridGetCoord", error)

   if (localpet==0) print*,"- CALL GridGetCoord FOR INPUT U GRID Y-COORD."
   nullify(lat_src_ptr)
   call ESMF_GridGetCoord(target_grid, &
                          staggerLoc=ESMF_STAGGERLOC_EDGE1, &
                          coordDim=2, &
                          computationalLBound=clb, &
                          computationalUBound=cub, &
                          farrayPtr=lat_src_ptr, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN GridGetCoord", error)

 ! Read coordinates from file and set on grid
 allocate(templat(clb(1):cub(1),clb(2):cub(2)))
 allocate(templon(clb(1):cub(1),clb(2):cub(2)))
 starts = (/clb(1),clb(2)/)
 counts = (/cub(1)-clb(1)+1,cub(2)-clb(2)+1/)
 
 if (localpet==0) print*,'- READ LONGITUDE U ID'
 error=nf90_inq_varid(ncid, 'XLONG_U', id_var)
 call netcdf_err(error, 'reading longitude u id')

 if (localpet==0) print*,'- READ LONGITUDE U'
 error=nf90_get_var(ncid, id_var, start=starts,count=counts,values=templon)
 call netcdf_err(error, 'reading longitude u')
 
 if (localpet==0) print*,'- READ LATITUDE U ID'
 error=nf90_inq_varid(ncid, 'XLAT_U', id_var)
 call netcdf_err(error, 'reading latitude u id')

 if (localpet==0) print*,'- READ LATITUDE U'
 error=nf90_get_var(ncid, id_var, start=starts,count=counts,values=templat)
 call netcdf_err(error, 'reading latitude u')
    
    do j = clb(2),cub(2)
      do i = clb(1), cub(1)
        lon_src_ptr(i,j)=real(templon(i,j),esmf_kind_r8)
        lat_src_ptr(i,j)=real(templat(i,j),esmf_kind_r8)
      enddo
    enddo

  nullify(lon_src_ptr)
  nullify(lat_src_ptr)

if (localpet==0) print*,"- CALL FieldCreate FOR TARGET GRID LATITUDE."
 latitude_u_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE1, &
                                   name="target_grid_latitude u", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", error)


 if (localpet==0) print*,"- CALL FieldCreate FOR TARGET GRID LONGITUDE."
 longitude_u_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE1, &
                                   name="target_grid_longitude u", &
                                   rc=error)

 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", error)

  call ESMF_FieldGet(latitude_u_target_grid, farrayPtr=lat_src_ptr,rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldGet", error)

  call ESMF_FieldGet(longitude_u_target_grid, farrayPtr=lon_src_ptr,rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldGet", error)

  do j = clb(2),cub(2)
      do i = clb(1), cub(1)
        lon_src_ptr(i,j)=real(templon(i,j),esmf_kind_r8)
        lat_src_ptr(i,j)=real(templat(i,j),esmf_kind_r8)
      enddo
    enddo
  nullify(lat_src_ptr)
  nullify(lon_src_ptr)
  deallocate(templat,templon)
  
!---------- N-S stagger grid coordinate creation

  if (localpet==0) print*,"- CALL GridGetCoord FOR INPUT V GRID X-COORD."
   nullify(lon_src_ptr)
   call ESMF_GridGetCoord(target_grid, &
                          staggerLoc=ESMF_STAGGERLOC_EDGE2, &
                          coordDim=1, &
                          farrayPtr=lon_src_ptr, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN GridGetCoord", error)

   if (localpet==0) print*,"- CALL GridGetCoord FOR INPUT V GRID Y-COORD."
   nullify(lat_src_ptr)
   call ESMF_GridGetCoord(target_grid, &
                          staggerLoc=ESMF_STAGGERLOC_EDGE2, &
                          coordDim=2, &
                          computationalLBound=clb, &
                          computationalUBound=cub, &
                          farrayPtr=lat_src_ptr, rc=error)
   if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN GridGetCoord", error)

 allocate(templat(clb(1):cub(1),clb(2):cub(2)))
 allocate(templon(clb(1):cub(1),clb(2):cub(2)))
 starts = (/clb(1),clb(2)/)
 counts = (/cub(1)-clb(1)+1,cub(2)-clb(2)+1/)
 

 if (localpet==0) print*,'- READ LONGITUDE ID'
 error=nf90_inq_varid(ncid, 'XLONG_V', id_var) ! CSS bug fix (was XLONGV)
 call netcdf_err(error, 'reading longitude id')


 if (localpet==0) print*,'- READ LONGITUDE V'
 error=nf90_get_var(ncid, id_var, start=starts,count=counts,values=templon)
 call netcdf_err(error, 'reading longitude v')
 
 if (localpet==0) print*,'- READ LATITUDE V ID'
 error=nf90_inq_varid(ncid, 'XLAT_V', id_var)
 call netcdf_err(error, 'reading latitude v id')

 if (localpet==0) print*,'- READ LATITUDE V'
 error=nf90_get_var(ncid, id_var, start=starts,count=counts,values=templat)
 call netcdf_err(error, 'reading latitude v')
    
    do j = clb(2),cub(2)
      do i = clb(1), cub(1)
        lon_src_ptr(i,j)=real(templon(i,j),esmf_kind_r8)
        lat_src_ptr(i,j)=real(templat(i,j),esmf_kind_r8)
      enddo
    enddo

  nullify(lon_src_ptr)
  nullify(lat_src_ptr)

if (localpet==0) print*,"- CALL FieldCreate FOR TARGET GRID LATITUDE."
 latitude_v_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE2, &
                                   name="target_grid_latitude_v", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", error)


 if (localpet==0) print*,"- CALL FieldCreate FOR TARGET GRID LONGITUDE."
 longitude_v_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE2, &
                                   name="target_grid_longitude_v", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", error)

  call ESMF_FieldGet(latitude_v_target_grid, farrayPtr=lat_src_ptr,rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldGet", error)

  call ESMF_FieldGet(longitude_v_target_grid, farrayPtr=lon_src_ptr,rc=error)
     if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldGet", error)

  do j = clb(2),cub(2)
      do i = clb(1), cub(1)
        lon_src_ptr(i,j)=real(templon(i,j),esmf_kind_r8)
        lat_src_ptr(i,j)=real(templat(i,j),esmf_kind_r8)
      enddo
    enddo
  nullify(lat_src_ptr)
  nullify(lon_src_ptr)

  deallocate(templat,templon)
  
!----------------------------------------------------------
!  Read in map factors and set on grid 
!----------------------------------------------------------

  if (localpet==0) print*,"- CALL FieldCreate FOR TARGET GRID mapfac_m."
 mapfac_m_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="target_grid_hgt", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", error)

 if (localpet==0) print*,"- CALL FieldGet FOR TARGET GRID mapfac_m. "
 call ESMF_FieldGet(mapfac_m_target_grid, farrayPtr=mapptr, &
                    computationalLBound=clb,computationalUBound=cub, &
                    rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
   call error_handler("IN FieldGet", error)
   
 allocate(mapfac_temp(clb(1):cub(1),clb(2):cub(2)))
 starts = (/clb(1),clb(2)/)
 counts = (/cub(1)-clb(1)+1,cub(2)-clb(2)+1/)
 
 if (localpet==0) print*,'- READ MAPFAC_M ID'
 error=nf90_inq_varid(ncid, 'MAPFAC_M', id_var)
 call netcdf_err(error, 'reading MAPFAC_M id')

 if (localpet==0) print*,'- READ MAPFAC_M'
 error=nf90_get_var(ncid, id_var, start=starts,count=counts,values=mapfac_temp)
 call netcdf_err(error, 'reading MAPFAC_M')

 do j=clb(2),cub(2)
 do i=clb(1),cub(1)
        mapptr(i,j) = mapfac_temp(i,j)
 enddo
 enddo
 nullify(mapptr)
 deallocate(mapfac_temp)

!----- MAPFAC U
 if (localpet==0) print*,"- CALL FieldCreate FOR TARGET GRID mapfac_u."
 mapfac_u_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE1, &
                                   name="target_grid_hgt", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", error)

 if (localpet==0) print*,"- CALL FieldGet FOR TARGET GRID mapfac_u. "
 call ESMF_FieldGet(mapfac_u_target_grid, farrayPtr=mapptr, &
                    computationalLBound=clb,computationalUBound=cub, &
                     rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
   call error_handler("IN FieldGet", error)

  allocate(mapfac_temp(clb(1):cub(1),clb(2):cub(2)))
 starts = (/clb(1),clb(2)/)
 counts = (/cub(1)-clb(1)+1,cub(2)-clb(2)+1/)
 
 if (localpet==0) print*,'- READ MAPFAC_U ID'
 error=nf90_inq_varid(ncid, 'MAPFAC_U', id_var)
 call netcdf_err(error, 'reading MAPFAC_U id')

 if (localpet==0) print*,'- READ MAPFACU'
 error=nf90_get_var(ncid, id_var, start=starts,count=counts,values=mapfac_temp)
 call netcdf_err(error, 'reading MAPFAC_U')

 do j=clb(2),cub(2)
 do i=clb(1),cub(1)
        mapptr(i,j) = mapfac_temp(i,j)
 enddo
 enddo
 nullify(mapptr)
 deallocate(mapfac_temp)

!------ MAPFAC V
 if (localpet==0) print*,"- CALL FieldCreate FOR TARGET GRID mapfac_v."
 mapfac_v_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_EDGE2, &
                                   name="target_grid_hgt", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", error)

 if (localpet==0) print*,"- CALL FieldGet FOR TARGET GRID mapfac_v. "
 call ESMF_FieldGet(mapfac_v_target_grid, farrayPtr=mapptr, &
                    computationalLBound=clb,computationalUBound=cub, &
                    rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
   call error_handler("IN FieldGet", error)

 allocate(mapfac_temp(clb(1):cub(1),clb(2):cub(2)))
 starts = (/clb(1),clb(2)/)
 counts = (/cub(1)-clb(1)+1,cub(2)-clb(2)+1/)
 
 if (localpet==0) print*,'- READ MAPFAC_V ID'
 error=nf90_inq_varid(ncid, 'MAPFAC_V', id_var)
 call netcdf_err(error, 'reading MAPFAC_V id')

 if (localpet==0) print*,'- READ MAPFAC_V'
 error=nf90_get_var(ncid, id_var, start=starts,count=counts,values=mapfac_temp)
 call netcdf_err(error, 'reading MAPFAC_V')

 do j=clb(2),cub(2)
 do i=clb(1),cub(1)
        mapptr(i,j) = mapfac_temp(i,j)
 enddo
 enddo
 nullify(mapptr)
 deallocate(mapfac_temp)

 if (proj_code==PROJ_LC) then
   !--- Read in grid angle arrays and set on grid
    if (localpet==0) print*, "- CALL FieldCreate for TARGET GRID sinalpha."
    sina_target_grid = ESMF_FieldCreate(target_grid, &
                               typekind=ESMF_TYPEKIND_R8, &
                               staggerloc=ESMF_STAGGERLOC_CENTER, &
                               name="target_grid_sina", &
                               rc=error)

    if (localpet==0) print*,"- CALL FieldGet FOR TARGET GRID sina. "
    call ESMF_FieldGet(sina_target_grid, farrayPtr=mapptr, &
                    computationalLBound=clb,computationalUBound=cub, &
                    rc=error)
    if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldGet", error)

    allocate(dum2d(clb(1):cub(1),clb(2):cub(2)))
    starts = (/clb(1),clb(2)/)
    counts = (/cub(1)-clb(1)+1,cub(2)-clb(2)+1/)
    if (localpet==0) print*,'- READ SINALPHA ID'
    error=nf90_inq_varid(ncid, 'SINALPHA', id_var)
    call netcdf_err(error, 'reading SINALPHA id')

    if (localpet==0) print*,'- READ SINALPHA'
    error=nf90_get_var(ncid, id_var, start=starts,count=counts,values=dum2d)
    call netcdf_err(error, 'reading SINALPHA')

    do j=clb(2),cub(2)
    do i=clb(1),cub(1)
        mapptr(i,j) = dum2d(i,j)
    enddo
    enddo
    nullify(mapptr)
    deallocate(dum2d)

    if (localpet==0) print*, "- CALL FieldCreate for TARGET GRID cosalpha."
    cosa_target_grid = ESMF_FieldCreate(target_grid, &
                               typekind=ESMF_TYPEKIND_R8, &
                               staggerloc=ESMF_STAGGERLOC_CENTER, &
                               name="target_grid_cosa", &
                               rc=error)

    if (localpet==0) print*,"- CALL FieldGet FOR TARGET GRID sina. "
    call ESMF_FieldGet(cosa_target_grid, farrayPtr=mapptr, &
                    computationalLBound=clb,computationalUBound=cub, &
                    rc=error)
    if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
      call error_handler("IN FieldGet", error)

    allocate(dum2d(clb(1):cub(1),clb(2):cub(2)))
    starts = (/clb(1),clb(2)/)
    counts = (/cub(1)-clb(1)+1,cub(2)-clb(2)+1/)
    if (localpet==0) print*,'- READ COSALPHA ID'
    error=nf90_inq_varid(ncid, 'COSALPHA', id_var)
    call netcdf_err(error, 'reading COSALPHA id')

    if (localpet==0) print*,'- READ COSALPHA'
    error=nf90_get_var(ncid, id_var, start=starts,count=counts,values=dum2d)
    call netcdf_err(error, 'reading COSALPHA')

    do j=clb(2),cub(2)
    do i=clb(1),cub(1)
        mapptr(i,j) = dum2d(i,j)
    enddo
    enddo

    nullify(mapptr)
    deallocate(dum2d)
 endif

!------- Height field
! CSS changed from ESMF_STAGGERLOC_EDGE2 to ESMF_STAGGERLOC_CENTER
  if (localpet==0) print*,"- CALL FieldCreate FOR TARGET GRID hgt."
 hgt_target_grid = ESMF_FieldCreate(target_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name="target_grid_hgt", &
                                   rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
    call error_handler("IN FieldCreate", error)

 if (localpet==0) print*,"- CALL FieldGet FOR TARGET GRID hgt. "
 call ESMF_FieldGet(hgt_target_grid, farrayPtr=hgtptr, &
                    computationalLBound=clb,computationalUBound=cub, &
                    rc=error)
 if(ESMF_logFoundError(rcToCheck=error,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__))&
   call error_handler("IN FieldGet", error)

 allocate(dum2d(clb(1):cub(1),clb(2):cub(2)))
 starts = (/clb(1),clb(2)/)
 counts = (/cub(1)-clb(1)+1,cub(2)-clb(2)+1/)

 if (localpet==0) print*,'- READ HGT ID'
 error=nf90_inq_varid(ncid, 'HGT', id_var)
 if (error .ne. 0) then
   error=nf90_inq_varid(ncid, 'HGT_M', id_var)
   call netcdf_err(error, 'reading hgt id')
 endif

 if (localpet==0) print*,'- READ HGT'
 error=nf90_get_var(ncid, id_var, start=starts,count=counts, values=dum2d)
 call netcdf_err(error, 'reading hgt')
 
 do j=clb(2),cub(2)
 do i=clb(1),cub(1)
        hgtptr(i,j) = dum2d(i,j)
 enddo
 enddo
 nullify(hgtptr)
 deallocate(dum2d)

 error = nf90_close(ncid)

 end subroutine define_target_grid_file

 !> For grids with equal cell sizes (e.g., lambert conformal), get lat and on of the grid
!! cell corners
!!
!! @param [in]  latitude  2d array of grid latitude
!! @param [in]  longitude 2d array of grid longitude
!! @param [in]  dx grid cell size in meters
!! @param [in]  clb  computational lower bound as returned from ESMF_GridGetCoord
!! @param [in]  cub  computational upper bound as returned from ESMF_GridGetCoord
!! @param [out] latitude_sw  returned latitude at sw corner of grid cells
!! @param [out]elongitude_sw  returned longitude at sw corner of grid cells
!! @author Larissa Reames CIWRO/NSSL/FRDD

  subroutine get_cell_corners( latitude, longitude, latitude_sw, longitude_sw, dx,clb,cub)
  implicit none

  real(esmf_kind_r8), intent(in), pointer :: latitude(:,:)
  real(esmf_kind_r8), intent(inout), pointer   :: latitude_sw(:,:)
  real(esmf_kind_r8), intent(in),pointer    :: longitude(:,:)
  real(esmf_kind_r8), intent(inout), pointer   :: longitude_sw(:,:)
  real(esmf_kind_r8), intent(in)    :: dx !grid cell side size (m)

  integer, intent(in) :: clb(2), cub(2)

  real(esmf_kind_r8)                :: lat1, lon1, lat2, lon2, d, brng


  real(esmf_kind_r8), parameter    :: pi = 3.14159265359
  real(esmf_kind_r8), parameter    :: R =  6370000.0
  real(esmf_kind_r8), parameter    :: bearingInDegrees = 135.0

  integer                           :: i, j, e

  d = sqrt((dx**2.0_esmf_kind_r8)/2.0_esmf_kind_r8)

  do j = clb(2),cub(2)
   do i = clb(1), cub(1)
         if (j == jp1_target .and. i == ip1_target) then
             lat1 = latitude(i_target,j_target)  * ( pi / 180.0_esmf_kind_r8 )
             lon1 = longitude(i_target,j_target) * ( pi / 180.0_esmf_kind_r8 )
             brng = 315.0_esmf_kind_r8 * pi / 180.0_esmf_kind_r8
             lat2 = asin( sin( lat1 ) * cos( d / R ) + cos( lat1 ) * sin( d / R ) * cos( brng ) );
             lon2= lon1 + atan2( sin( brng ) * sin( d / R ) * cos( lat1 ), cos( d / R ) - sin( lat1 ) * sin( lat2 ) );
             latitude_sw(ip1_target,jp1_target) = lat2 * 180.0_esmf_kind_r8 / pi
             longitude_sw(ip1_target,jp1_target) = lon2 * 180.0_esmf_kind_r8 / pi
             cycle
         endif

     if (i == ip1_target) then
       brng = 225.0_esmf_kind_r8 * pi / 180.0_esmf_kind_r8
       lat1 = latitude(i_target,j)  * ( pi / 180.0_esmf_kind_r8 )
       lon1 = longitude(i_target,j) * ( pi / 180.0_esmf_kind_r8 )
       lat2 = asin( sin( lat1 ) * cos( d / R ) + cos( lat1 ) * sin( d / R ) * cos( brng ) );
       lon2= lon1 + atan2( sin( brng ) * sin( d / R ) * cos( lat1 ), cos( d / R ) - sin( lat1 ) * sin( lat2 ) );
       latitude_sw(ip1_target,j) = lat2 * 180.0_esmf_kind_r8 / pi
       longitude_sw(ip1_target,j) = lon2 * 180.0_esmf_kind_r8 / pi
       cycle
     endif

     if (j == jp1_target) then
       brng = 45.0_esmf_kind_r8 * pi / 180.0_esmf_kind_r8
       lat1 = latitude(i,j_target)  * ( pi / 180.0_esmf_kind_r8 )
       lon1 = longitude(i,j_target) * ( pi / 180.0_esmf_kind_r8 )
       lat2 = asin( sin( lat1 ) * cos( d / R ) + cos( lat1 ) * sin( d / R ) * cos( brng ) );
       lon2= lon1 + atan2( sin( brng ) * sin( d / R ) * cos( lat1 ), cos( d / R ) - sin( lat1 ) * sin( lat2 ) );
       latitude_sw(i,jp1_target) = lat2 * 180.0_esmf_kind_r8 / pi
       longitude_sw(i,jp1_target) = lon2 * 180.0_esmf_kind_r8 / pi
       cycle
     endif

     lat1 = latitude(i,j)  * ( pi / 180.0_esmf_kind_r8 )
     lon1 = longitude(i,j) * ( pi / 180.0_esmf_kind_r8 )

     brng = bearingInDegrees * ( pi / 180.0_esmf_kind_r8 );
     lat2 = asin( sin( lat1 ) * cos( d / R ) + cos( lat1 ) * sin( d / R ) * cos( brng ) );
     lon2= lon1 + atan2( sin( brng ) * sin( d / R ) * cos( lat1 ), cos( d / R ) - sin( lat1 ) * sin( lat2 ) );

     latitude_sw(i,j) = lat2 * 180.0_esmf_kind_r8 / pi
     longitude_sw(i,j) = lon2 * 180.0_esmf_kind_r8 / pi

   enddo
 enddo

 end subroutine get_cell_corners

 subroutine cleanup_input_target_grid_data(localpet)

 use program_setup, only    : interp_diag, interp_hist


 implicit none
 integer, intent(in)              :: localpet
 integer                          :: rc, i
 type(esmf_field), allocatable    :: fields(:)

 if (localpet==0) print*,"- DESTROY MODEL DATA."

 call ESMF_FieldDestroy(node_latitude_input_grid,rc=rc)
 call ESMF_FieldDestroy(node_longitude_input_grid,rc=rc)

 call ESMF_FieldDestroy(cell_latitude_input_grid,rc=rc)
 call ESMF_FieldDestroy(cell_longitude_input_grid,rc=rc)

 call ESMF_FieldDestroy(latitude_target_grid, rc=rc)
 call ESMF_FieldDestroy(longitude_target_grid, rc=rc)

 call ESMF_FieldDestroy(latitude_u_target_grid, rc=rc)
 call ESMF_FieldDestroy(longitude_u_target_grid, rc=rc)

 call ESMF_FieldDestroy(latitude_v_target_grid, rc=rc)
 call ESMF_FieldDestroy(longitude_v_target_grid, rc=rc)

 if (interp_diag) then
    allocate(fields(n_diag_fields))
    call ESMF_FieldBundleGet(input_diag_bundle, fieldList=fields, &
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=rc)
    do i = 1, n_diag_fields
        call ESMF_FieldDestroy(fields(i), rc=rc)
    enddo

    call ESMF_FieldBundleDestroy(input_diag_bundle)

    call ESMF_FieldBundleGet(target_diag_bundle, fieldList=fields, &
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=rc)
    do i = 1, n_diag_fields
        call ESMF_FieldDestroy(fields(i), rc=rc)
    enddo

    call ESMF_FieldBundleDestroy(target_diag_bundle)
    deallocate(fields)
 endif

 if (n_hist_fields_2d_cons>0) then
    allocate(fields(n_hist_fields_2d_cons))
    call ESMF_FieldBundleGet(input_hist_bundle_2d_cons, fieldList=fields, &
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=rc)
    do i = 1, n_hist_fields_2d_cons
        call ESMF_FieldDestroy(fields(i), rc=rc)
    enddo

    call ESMF_FieldBundleDestroy(input_hist_bundle_2d_cons)

    call ESMF_FieldBundleGet(target_hist_bundle_2d_cons, fieldList=fields, &
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=rc)
    do i = 1, n_hist_fields_2d_cons
        call ESMF_FieldDestroy(fields(i), rc=rc)
    enddo

    call ESMF_FieldBundleDestroy(target_hist_bundle_2d_cons)
    deallocate(fields)
 endif

  if (n_hist_fields_2d_nstd>0) then
    allocate(fields(n_hist_fields_2d_nstd))
    call ESMF_FieldBundleGet(input_hist_bundle_2d_nstd, fieldList=fields, &
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=rc)
    do i = 1, n_hist_fields_2d_nstd
        call ESMF_FieldDestroy(fields(i), rc=rc)
    enddo

    call ESMF_FieldBundleDestroy(input_hist_bundle_2d_nstd)

    call ESMF_FieldBundleGet(target_hist_bundle_2d_nstd, fieldList=fields, &
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=rc)
    do i = 1, n_hist_fields_2d_nstd
        call ESMF_FieldDestroy(fields(i), rc=rc)
    enddo

    call ESMF_FieldBundleDestroy(target_hist_bundle_2d_nstd)
    deallocate(fields)
 endif

 if (n_hist_fields_2d_patch>0) then
    allocate(fields(n_hist_fields_2d_patch))
    call ESMF_FieldBundleGet(input_hist_bundle_2d_patch, fieldList=fields, &
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=rc)
    do i = 1, n_hist_fields_2d_patch
        call ESMF_FieldDestroy(fields(i), rc=rc)
    enddo

    call ESMF_FieldBundleDestroy(input_hist_bundle_2d_patch)

    call ESMF_FieldBundleGet(target_hist_bundle_2d_patch, fieldList=fields, &
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=rc)
    do i = 1, n_hist_fields_2d_patch
        call ESMF_FieldDestroy(fields(i), rc=rc)
    enddo

    call ESMF_FieldBundleDestroy(target_hist_bundle_2d_patch)
    deallocate(fields)
 endif

  if (n_hist_fields_3d_nz>0) then
    allocate(fields(n_hist_fields_3d_nz))
    call ESMF_FieldBundleGet(input_hist_bundle_3d_nz, fieldList=fields, &
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=rc)
    do i = 1, n_hist_fields_3d_nz
        call ESMF_FieldDestroy(fields(i), rc=rc)
    enddo

    call ESMF_FieldBundleDestroy(input_hist_bundle_3d_nz)

    call ESMF_FieldBundleGet(target_hist_bundle_3d_nz, fieldList=fields, &
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=rc)
    do i = 1, n_hist_fields_3d_nz
        call ESMF_FieldDestroy(fields(i), rc=rc)
    enddo

    call ESMF_FieldBundleDestroy(target_hist_bundle_3d_nz)
    deallocate(fields)
 endif

 if (n_hist_fields_3d_nzp1>0) then
    allocate(fields(n_hist_fields_3d_nzp1))
    call ESMF_FieldBundleGet(input_hist_bundle_3d_nzp1, fieldList=fields, &
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=rc)
    do i = 1, n_hist_fields_3d_nzp1
        call ESMF_FieldDestroy(fields(i), rc=rc)
    enddo

    call ESMF_FieldBundleDestroy(input_hist_bundle_3d_nzp1)

    call ESMF_FieldBundleGet(target_hist_bundle_3d_nzp1, fieldList=fields, &
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=rc)
    do i = 1, n_hist_fields_3d_nzp1
        call ESMF_FieldDestroy(fields(i), rc=rc)
    enddo

    call ESMF_FieldBundleDestroy(target_hist_bundle_3d_nzp1)
    deallocate(fields)
 endif

 if (n_hist_fields_3d_vert>0) then
    allocate(fields(n_hist_fields_3d_vert))
    call ESMF_FieldBundleGet(input_hist_bundle_3d_vert, fieldList=fields, &
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=rc)
    do i = 1, n_hist_fields_3d_vert
        call ESMF_FieldDestroy(fields(i), rc=rc)
    enddo

    call ESMF_FieldBundleDestroy(input_hist_bundle_3d_vert)

    call ESMF_FieldBundleGet(target_hist_bundle_3d_vert, fieldList=fields, &
                          itemorderflag=ESMF_ITEMORDER_ADDORDER, &
                          rc=rc)
    do i = 1, n_hist_fields_3d_vert
        call ESMF_FieldDestroy(fields(i), rc=rc)
    enddo

    call ESMF_FieldBundleDestroy(target_hist_bundle_3d_vert)
    deallocate(fields)
 endif
 call ESMF_MeshDestroy(input_grid, rc=rc)

 call ESMF_GridDestroy(target_grid, rc=rc)

 end subroutine cleanup_input_target_grid_data

subroutine unique_sort(val,nvals, final)
    implicit none
    integer, intent(in) :: nvals
    integer, intent(in) :: val(nvals)
    integer, intent(inout), allocatable :: final(:)
    integer :: i = 0, min_val, max_val
    integer, dimension(nvals) :: unique

    min_val = minval(val)-1
    max_val = maxval(val)
    do while (min_val<max_val)
        min_val = minval(val, mask=val>min_val)
        if (min_val>0) then
          i = i+1
          unique(i) = min_val
        endif
    enddo
    allocate(final(i), source=unique(1:i))   !<-- Or, just use unique(1:i)
end subroutine unique_sort

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: get_lat_lon_fields
   !
   ! Purpose: To calculate the latitude and longitude for every gridpoint in the
   !   tile of the model domain. The caller may specify that the grid for which 
   !   values are computed is staggered or unstaggered using the "stagger" 
   !   argument.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_lat_lon_fields(xlat_arr, xlon_arr, start_mem_i, &
                                 start_mem_j, end_mem_i, end_mem_j, stagger, comp_ll, &
                                 sub_x, sub_y)
   
      use llxy_module
    
      implicit none
    
      ! Arguments
      integer, intent(in) :: start_mem_i, start_mem_j, end_mem_i, &
                             end_mem_j, stagger
      real, dimension(start_mem_i:end_mem_i, start_mem_j:end_mem_j), intent(out) :: xlat_arr, xlon_arr
      logical, optional, intent(in) :: comp_ll
      integer, optional, intent(in) :: sub_x, sub_y

      ! Local variables
      integer :: i, j
      real :: rx, ry
    
      rx = 1.0
      ry = 1.0
      if (present(sub_x)) rx = real(sub_x)
      if (present(sub_y)) ry = real(sub_y)

      do i=start_mem_i, end_mem_i
         do j=start_mem_j, end_mem_j
            call xytoll(real(i-0.5)/rx+0.5, real(j-0.5)/ry+0.5, &
                        xlat_arr(i,j), xlon_arr(i,j), stagger, comp_ll=comp_ll)
         end do
      end do

   end subroutine get_lat_lon_fields

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: get_map_factor
   !
   ! Purpose: Given the latitude field, this routine calculates map factors for 
   !   the grid points of the specified domain. For different grids (e.g., C grid, 
   !   E grid), the latitude array should provide the latitudes of the points for
   !   which map factors are to be calculated. 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_map_factor(xlat_arr, xlon_arr, mapfac_arr_x, mapfac_arr_y, &
                             start_mem_i, start_mem_j, end_mem_i, end_mem_j)
   
      implicit none
    
      ! Arguments
      integer, intent(in) :: start_mem_i, start_mem_j, end_mem_i, end_mem_j
      real, dimension(start_mem_i:end_mem_i, start_mem_j:end_mem_j), intent(in) :: xlat_arr, xlon_arr
      real, dimension(start_mem_i:end_mem_i, start_mem_j:end_mem_j), intent(out) :: mapfac_arr_x
      real, dimension(start_mem_i:end_mem_i, start_mem_j:end_mem_j), intent(out) :: mapfac_arr_y
    
      ! Local variables
      integer :: i, j
      real :: n, colat, colat0, colat1, colat2, comp_lat, comp_lon
    
      !
      ! Equations for map factor given in Principles of Meteorological Analysis,
      ! Walter J. Saucier, pp. 32-33 
      !
    
      ! Lambert conformal projection
      if (proj_code == PROJ_LC) then
         if (truelat1 /= truelat2) then
            colat1 = rad_per_deg*(90.0 - truelat1)
            colat2 = rad_per_deg*(90.0 - truelat2)
            n = (log(sin(colat1)) - log(sin(colat2))) &
                / (log(tan(colat1/2.0)) - log(tan(colat2/2.0)))
      
            do i=start_mem_i, end_mem_i
               do j=start_mem_j, end_mem_j
                  colat = rad_per_deg*(90.0 - xlat_arr(i,j))
                  mapfac_arr_x(i,j) = sin(colat2)/sin(colat)*(tan(colat/2.0)/tan(colat2/2.0))**n
                  mapfac_arr_y(i,j) = mapfac_arr_x(i,j)
               end do
            end do
     
         else
            colat0 = rad_per_deg*(90.0 - truelat1)
      
            do i=start_mem_i, end_mem_i
               do j=start_mem_j, end_mem_j
                  colat = rad_per_deg*(90.0 - xlat_arr(i,j))
                  mapfac_arr_x(i,j) = sin(colat0)/sin(colat)*(tan(colat/2.0)/tan(colat0/2.0))**cos(colat0)
                  mapfac_arr_y(i,j) = mapfac_arr_x(i,j)
               end do
            end do
    
         end if
    
      ! Polar stereographic projection
      else if (proj_code == PROJ_PS) then
    
         do i=start_mem_i, end_mem_i
            do j=start_mem_j, end_mem_j
               mapfac_arr_x(i,j) = (1.0 + sin(rad_per_deg*abs(truelat1)))/(1.0 + sin(rad_per_deg*sign(1.,truelat1)*xlat_arr(i,j)))
               mapfac_arr_y(i,j) = mapfac_arr_x(i,j)
            end do
         end do
    
      ! Mercator projection 
      else if (proj_code == PROJ_MERC) then
         colat0 = rad_per_deg*(90.0 - truelat1)
     
         do i=start_mem_i, end_mem_i
            do j=start_mem_j, end_mem_j
               colat = rad_per_deg*(90.0 - xlat_arr(i,j))
               mapfac_arr_x(i,j) = sin(colat0) / sin(colat) 
               mapfac_arr_y(i,j) = mapfac_arr_x(i,j)
            end do
         end do
    
      ! Global cylindrical projection
      else if (proj_code == PROJ_CYL) then
     
         do i=start_mem_i, end_mem_i
            do j=start_mem_j, end_mem_j
               if (abs(xlat_arr(i,j)) == 90.0) then
                  mapfac_arr_x(i,j) = 0.    ! MSF actually becomes infinite at poles, but 
                                            !   the values should never be used there; by
                                            !   setting to 0, we hope to induce a "divide
                                            !   by zero" error if they are
               else
                  mapfac_arr_x(i,j) = 1.0 / cos(xlat_arr(i,j)*rad_per_deg) 
               end if
               mapfac_arr_y(i,j) = 1.0
            end do
         end do
    
      ! Rotated global cylindrical projection
      else if (proj_code == PROJ_CASSINI) then
     
         if (abs(pole_lat) == 90.) then
            do i=start_mem_i, end_mem_i
               do j=start_mem_j, end_mem_j
                  if (abs(xlat_arr(i,j)) >= 90.0) then
                     mapfac_arr_x(i,j) = 0.    ! MSF actually becomes infinite at poles, but 
                                               !   the values should never be used there; by
                                               !   setting to 0, we hope to induce a "divide
                                               !   by zero" error if they are
                  else
                     mapfac_arr_x(i,j) = 1.0 / cos(xlat_arr(i,j)*rad_per_deg) 
                  end if
                  mapfac_arr_y(i,j) = 1.0
               end do
            end do
         else
            do i=start_mem_i, end_mem_i
               do j=start_mem_j, end_mem_j
                  call rotate_coords(xlat_arr(i,j),xlon_arr(i,j), &
                                     comp_lat, comp_lon, &
                                     pole_lat, pole_lon, stand_lon, &
                                     -1)
                  if (abs(comp_lat) >= 90.0) then
                     mapfac_arr_x(i,j) = 0.    ! MSF actually becomes infinite at poles, but 
                                               !   the values should never be used there; by
                                               !   setting to 0, we hope to induce a "divide
                                               !   by zero" error if they are
                  else
                     mapfac_arr_x(i,j) = 1.0 / cos(comp_lat*rad_per_deg) 
                  end if
                  mapfac_arr_y(i,j) = 1.0
               end do
            end do
         end if
    
      else if (proj_code == PROJ_ROTLL) then
    
         do i=start_mem_i, end_mem_i
            do j=start_mem_j, end_mem_j
               mapfac_arr_x(i,j) = 1.0
               mapfac_arr_y(i,j) = 1.0
            end do
         end do
    
      end if
   
   end subroutine get_map_factor

   subroutine read_block_decomp_file(localpet, npets,file,ncells,myCells,myCells_num)

    implicit none

   integer, INTENT(IN)                      :: localpet, ncells, npets
   character(500), INTENT(IN)              :: file
   integer, INTENT(OUT), ALLOCATABLE        :: myCells(:)
   integer, INTENT(OUT)                     :: myCells_num
   logical                                  :: file_exists
   integer :: k, istat, nlines, proc, proc_max
   integer, ALLOCATABLE                     :: myCells_temp(:)
   character(200) :: line,msg
   character(100) :: str1,str2

   INQUIRE(FILE=file, EXIST=file_exists)
   
  if (.not. file_exists) then
     call error_handler("BLOCK DECOMP FILE DOES NOT EXIST",-1)
   endif

   open(14, file=trim(file), form='formatted', iostat=istat)
   if (istat /= 0) then
     call error_handler("OPENING BLOCK DECOMP FILE", istat)
   endif

   nlines = 0

   !Loop over lines of file to count the number of cells
   do
     read(14, '(A)', iostat=istat) line
     if (istat/=0) exit
     if ( trim(line) .eq. '' ) cycle
     nlines = nlines+1
   enddo
   if (nlines /= ncells) call error_handler("BLOCK DECOMPOSITION FILE CONTAINS MORE CELLS THAN INPUT GRID", -1)

   allocate(myCells_temp(ncells))
  
   myCells_num = 0
   proc_max=0
   rewind(14)
    do k = 1,nlines
      read(14, *, iostat=istat) proc
     if (istat /= 0) call error_handler("READING BLOCK DECOMPOSITION FILE", istat)
     proc_max = max(proc,proc_max)
     if (localpet==proc) then
        myCells_num = myCells_num+1 
        myCells_temp(myCells_num) = k
     endif
    enddo
    
    write(str1,'(A,I10,A)') "BLOCK DECOMPOSITION FILE GENERATED FOR ",proc_max+1," PROCESSES BUT "
    write(str2,'(I10,A)') npets," PROCESSORS USED."
    write(msg,'(A)') trim(str1)//new_line('A')//trim(str2)
    if(proc_max+1 /= npets) call error_handler(trim(msg),-1)
    allocate(myCells(myCells_num))
    myCells(1:myCells_num) = myCells_temp(1:myCells_num)
   close(14)

   end subroutine read_block_decomp_file

   subroutine para_range(n1, n2, nprocs, irank, ista, iend)

      integer, intent(in) :: n1, n2, nprocs, irank
      integer, intent(out) :: ista, iend

      integer :: iwork1, iwork2

      iwork1 = (n2 - n1 + 1) / nprocs
      iwork2 = mod(n2 - n1 + 1, nprocs)
      ista = irank * iwork1 + n1 + min(irank, iwork2)
      iend = ista + iwork1 - 1
      if (iwork2 > irank) iend = iend + 1
      return
   end subroutine para_range
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: get_rotang
   !
   ! Purpose: To calculate the sine and cosine of rotation angle.
   !
   ! NOTES: The formulas used in this routine come from those in the
   !   vecrot_rotlat() routine of the original WRF SI.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_rotang(xlat_arr, xlon_arr, cosa, sina, &
                         start_mem_i, start_mem_j, end_mem_i, end_mem_j)

      implicit none

     ! Arguments
      integer, intent(in) :: start_mem_i, start_mem_j, end_mem_i, end_mem_j
      real(esmf_kind_r8), dimension(start_mem_i:end_mem_i, start_mem_j:end_mem_j), intent(in) :: xlat_arr, xlon_arr
      real, pointer, dimension(:,:), intent(inout) :: cosa, sina
      ! Local variables
      integer :: i, j
      real :: alpha, d_lon

      do i=start_mem_i, end_mem_i
         do j=start_mem_j+1, end_mem_j-1
            d_lon = xlon_arr(i,j+1)-xlon_arr(i,j-1)
            if (d_lon > 180.) then
               d_lon = d_lon - 360.
            else if (d_lon < -180.) then
               d_lon = d_lon + 360.
            end if

            alpha = atan2(-cos(xlat_arr(i,j)*RAD_PER_DEG) * (d_lon*RAD_PER_DEG), &
                            ((xlat_arr(i,j+1)-xlat_arr(i,j-1))*RAD_PER_DEG))
            sina(i,j) = sin(alpha)
            cosa(i,j) = cos(alpha)
         end do
      end do

      do i=start_mem_i, end_mem_i
         d_lon = xlon_arr(i,start_mem_j+1)-xlon_arr(i,start_mem_j)
         if (d_lon > 180.) then
            d_lon = d_lon - 360.
         else if (d_lon < -180.) then
            d_lon = d_lon + 360.
         end if

         alpha = atan2(-cos(xlat_arr(i,start_mem_j)*RAD_PER_DEG) * (d_lon*RAD_PER_DEG), &
                       ((xlat_arr(i,start_mem_j+1)-xlat_arr(i,start_mem_j))*RAD_PER_DEG))
         sina(i,start_mem_j) = sin(alpha)
         cosa(i,start_mem_j) = cos(alpha)
      end do

      do i=start_mem_i, end_mem_i
         d_lon = xlon_arr(i,end_mem_j)-xlon_arr(i,end_mem_j-1)
         if (d_lon > 180.) then
            d_lon = d_lon - 360.
         else if (d_lon < -180.) then
            d_lon = d_lon + 360.
         end if

         alpha = atan2(-cos(xlat_arr(i,end_mem_j)*RAD_PER_DEG) * (d_lon*RAD_PER_DEG), &
                       ((xlat_arr(i,end_mem_j)-xlat_arr(i,end_mem_j-1))*RAD_PER_DEG))
         sina(i,end_mem_j) = sin(alpha)
         cosa(i,end_mem_j) = cos(alpha)
      end do

   end subroutine get_rotang

 end module model_grid
