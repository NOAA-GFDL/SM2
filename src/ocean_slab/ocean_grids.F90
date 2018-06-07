module ocean_grids_mod
!
!<CONTACT EMAIL="Matthew.Harrison@noaa.gov"> M.J. Harrison 
!</CONTACT>
!
!<CONTACT EMAIL="Seth.Underwood@noaa.gov"> S. Underwood 
!</CONTACT>
!
!<CONTACT EMAIL="Zhi.Liang@noaa.gov"> Zhi Liang 
!</CONTACT>
!
! <REVIEWER EMAIL="">
! </REVIEWER>
!
!<OVERVIEW>
! Set up the ocean model grid spacing 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module sets up the ocean model grid based on information read in 
! from the grid_spec.nc file. It translates the generic names from the 
! grid_spec.nc file to the names used by the mixed layer ocean. 
!</DESCRIPTION>
!
!<NAMELIST NAME="ocean_grids_nml">
!  <DATA NAME="debug_grid" TYPE="logical">
!  For debugging. Note that most of the debugging stuff 
!  has been removed, but keep flag around in case need in future.
!  </DATA> 
!  <DATA NAME="verbose_init" TYPE="logical">
!  Prints out lots of initial checksums.  Useful to have on, so 
!  defaulted to true. 
!   </DATA> 
! </NAMELIST>
!
use fms_mod,           only: write_version_number
use fms_mod,           only: read_data, write_data
use fms_mod,           only: open_namelist_file, close_file, check_nml_error
use fms_mod,           only: field_exist, field_size, get_global_att_value
use mpp_domains_mod,   only: mpp_update_domains, mpp_global_field, domain2d
use mpp_domains_mod,   only: mpp_global_sum, BITWISE_EXACT_SUM
use mpp_mod,           only: stdout, stdlog, mpp_error, mpp_chksum, mpp_sum, FATAL, NOTE
use mpp_mod,           only: mpp_pe, mpp_root_pe
use mosaic_mod,        only: get_mosaic_ntiles, get_mosaic_ncontacts, get_mosaic_contact
use diag_output_mod,   only: get_diag_global_att, set_diag_global_att

implicit none

private

public ocean_grids_init
public set_ocean_grid_size
public set_ocean_hgrid_arrays


character(len=128), parameter :: version = '$Id$'
character(len=128), parameter :: tagname = '$Name$'

logical :: module_is_initialized = .FALSE.

!mosaic         
integer             :: grid_version 
integer, parameter  :: VERSION_0 = 0  ! grid file with field geolon_c
integer, parameter  :: VERSION_1 = 1  ! grid file with field x_C
integer, parameter  :: VERSION_2 = 2  ! mosaic file

! nml variables 
logical :: debug_this_module = .false.
logical :: verbose_init = .false.

! Grid variables
character(len=256) :: name
integer :: grid_ni, grid_nj
real, dimension(:,:), pointer :: x_c, y_c, grid_wet
real, dimension(:), pointer :: grid_x_c, grid_y_c
logical :: tripolar, mosaic

!type(ocean_grid_type) :: Grid

namelist /ocean_grids_nml/ debug_this_module, verbose_init

! grid file name
character(len=256) :: ocean_hgrid        ! will be set in set_ocean_grid_size
character(len=256) :: ocean_vgrid = 'INPUT/ocean_vgrid.nc'

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_grids_init">
!
! <DESCRIPTION>
! Initialize the grids module. 
! </DESCRIPTION>
!
subroutine ocean_grids_init(debug)

  logical, intent(in), optional :: debug
  integer                       :: ioun, io_status, ierr
  integer                       :: unit ! needed for writes to stdlog or stdout on DOE systems
  if ( module_is_initialized ) then
    call mpp_error(FATAL, '==>Error from ocean_grids_mod (ocean_grids_init): module has been initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  if (PRESENT(debug)) debug_this_module = debug

  ! provide for namelist over-ride of defaults 
  ioun = open_namelist_file()
  read (ioun,ocean_grids_nml,IOSTAT=io_status)
  if ( mpp_pe().EQ.mpp_root_pe() ) then 
     write (stdout(),'(/)')
     write (stdout(),ocean_grids_nml)
  end if
  unit = stdlog()
  write (unit,ocean_grids_nml)
  ierr = check_nml_error(io_status, 'ocean_grids_nml')
  call close_file(ioun)


end subroutine ocean_grids_init
! </SUBROUTINE> NAME="ocean_grids_init"


!#######################################################################
! <SUBROUTINE NAME="set_ocean_grid_size">
!
! <DESCRIPTION>
! Set the ocean grid size.  Model expects the grid specification file
! to be called grid_spec.nc.  
! </DESCRIPTION>
!
subroutine set_ocean_grid_size(ni,nj,grid_file, grid_name)

  character(len=*),      intent(in), optional :: grid_file
  character(len=*),      intent(in), optional :: grid_name
  integer, intent(inout)                      :: ni, nj
  
  integer                                   :: siz(4)
  integer                                   :: m
  integer                                   :: nx(1), ny(1)
  integer                                   :: ntiles, ncontacts
  integer, dimension(2)                     :: tile1, tile2
  integer, dimension(2)                     :: istart1, iend1, jstart1, jend1
  integer, dimension(2)                     :: istart2, iend2, jstart2, jend2
  character(len=256)                        :: grd_file, ocean_mosaic, attvalue

  
  grid_ni=0 ; grid_nj=0 
  mosaic = .false.
  grd_file = "INPUT/grid_spec.nc"
  if(present(grid_file)) grd_file = grid_file

  !
  !  Determine if the grid is mosaic file
  !
  if(field_exist(grd_file, 'ocn_mosaic_file') .or. field_exist(grd_file, 'gridfiles') ) then ! read from mosaic file
     if ( mpp_pe().EQ.mpp_root_pe() ) write(stdout(),*) '==>Note from ocean_grids_mod(set_ocean_grid_size): read grid from mosaic version grid'
     grid_version = VERSION_2
     mosaic = .true.
     if( field_exist(grd_file, 'ocn_mosaic_file') ) then ! coupler mosaic
        call read_data(grd_file, "ocn_mosaic_file", ocean_mosaic)
        ocean_mosaic = "INPUT/"//trim(ocean_mosaic)
     else
        ocean_mosaic = trim(grd_file)
     end if
     ntiles = get_mosaic_ntiles(ocean_mosaic)
     if(ntiles .NE. 1) call mpp_error(FATAL, '==>Error from ocean_grids_mod(set_ocean_grid_size): '//&
                   'ntiles should be 1 for ocean mosaic, contact developer')
! Notify diag_manager of mosaic grid
     call set_diag_global_att ('ocean','mosaic','1')
     call read_data(ocean_mosaic, "gridfiles", ocean_hgrid)
     ocean_hgrid = 'INPUT/'//trim(ocean_hgrid)
     call field_size(ocean_hgrid, 'x', siz)
     grid_ni = nx(1)
     grid_nj = ny(1)
     if(mod(siz(1),2) .NE. 1) call mpp_error(FATAL, '==>Error from ocean_grids_mod(set_ocean_grid_size): '//&
          'x-size of x in file '//trim(ocean_hgrid)//' should be 2*ni+1')
     if(mod(siz(2),2) .NE. 1) call mpp_error(FATAL, '==>Error from ocean_grids_mod(set_ocean_grid_size): '//&
          'y-size of x in file '//trim(ocean_hgrid)//' should be 2*nj+1')
     grid_ni = siz(1)/2 
     grid_nj = siz(2)/2 
  else  if(field_exist(grd_file, 'x_C')) then
     ocean_hgrid = grd_file
     if ( mpp_pe().EQ.mpp_root_pe()) write(stdout(),*) '==>Note from ocean_grids_mod(set_ocean_grid_size): read grid from new version grid'
     grid_version = VERSION_1
     call field_size( ocean_hgrid, 'x_C', siz)
     grid_ni = siz(1)
     grid_nj = siz(2) 
  else if(field_exist(grd_file, 'geolon_c')) then
     ocean_hgrid = grd_file
     if ( mpp_pe().EQ.mpp_root_pe() ) write(stdout(),*) '==>Note from ocean_grids_mod(set_ocean_grid_size): read grid from old version grid'
     grid_version = VERSION_0 
     call field_size( ocean_hgrid, 'geolon_c', siz)  
     grid_ni = siz(1)
     grid_nj = siz(2)
  else
     call mpp_error(FATAL, '==>Error from ocean_grids_mod(set_ocean_grid_size): '//&
                     'x_C, geolon_c, ocn_mosaic_file, gridfiles does not exist in file ' //trim(grd_file))
  endif

  if (grid_ni == 0 .or. grid_nj == 0 ) then
     if ( mpp_pe().EQ.mpp_root_pe() ) write(stdout(),*) '==>Error reading grid information from ',trim(grd_file),'. Make sure file exists'
     call mpp_error(FATAL,'==>Error reading grid information from grid file.  Are you sure file exists?')
  endif

  tripolar=.false.


  if(grid_version == VERSION_2) then
     !z1l: f_plane, beta_plane area not supported in mosaic grid. Need to think about to implement this. 
     if(field_exist(ocean_mosaic, "contacts") ) then
        ncontacts = get_mosaic_ncontacts(ocean_mosaic)
        if(ncontacts < 1) call mpp_error(FATAL,'==>Error from ocean_grids_mod(set_ocean_grid_size): '//&
                       'number of contacts should be larger than 0 when field contacts exist in file '//trim(ocean_mosaic) )
        if(ncontacts > 2) call mpp_error(FATAL,'==>Error from ocean_grids_mod(set_ocean_grid_size): '//&
                       'number of contacts should be no larger than 2')
        call get_mosaic_contact( ocean_mosaic, tile1(1:ncontacts), tile2(1:ncontacts),           &
             istart1(1:ncontacts), iend1(1:ncontacts), jstart1(1:ncontacts), jend1(1:ncontacts), &
             istart2(1:ncontacts), iend2(1:ncontacts), jstart2(1:ncontacts), jend2(1:ncontacts)  )
        do m = 1, ncontacts
           if(istart1(m) == iend1(m) ) then  ! x-direction contact, only cyclic condition
              if(istart2(m) .NE. iend2(m) ) call mpp_error(FATAL,  &
                   "==>Error from ocean_grids_mod(set_ocean_grid_size): only cyclic condition is allowed for x-boundary")
           else if( jstart1(m) == jend1(m) ) then  ! y-direction contact, cyclic or folded-north
              if(jstart2(m) .NE. jend2(m) ) call mpp_error(FATAL,  &
                   "==>Error from ocean_grids_mod(set_ocean_grid_size): "//&
                   "only cyclic/folded-north condition is allowed for y-boundary")
              if( jstart1(m) == jstart2(m) ) then ! folded north
             tripolar = .true.
         endif
           else 
              call mpp_error(FATAL,  &
                   "==>Error from ocean_grids_mod(set_ocean_grid_size): invalid boundary contact")
         end if
  end do
     end if
  else
     if( get_global_att_value(ocean_hgrid, "y_boundary_type", attvalue) ) then
        if(attvalue == 'fold_north_edge') then
           tripolar = .true.
        end if
     end if
  end if


  if(tripolar) then
     call mpp_error(NOTE,'==>Note from ocean_grids_mod(set_ocean_grid_size): y_boundary_type is fold_north_edge')
  endif



  if (PRESENT(grid_name)) then
     name = grid_name
  else
     name = 'ocean'
  endif

  ni = grid_ni
  nj = grid_nj
  
end subroutine set_ocean_grid_size
! </SUBROUTINE> NAME="set_ocean_grid_size"
!
!Following are the available fields in mosaic files
!--------------------------------------------------------
!Mosaic file     fields
!--------------------------------------------------------
!ocean_hgrid.nc  x, y, dx, dy, angle_dx, area
!ocean_vgrid.nc  zeta
!topog.nc        depth
!
! </DESCRIPTION>
subroutine set_ocean_hgrid_arrays(geolon,geolat,wet,domain)

  real, dimension(:,:), intent(inout) :: geolon, geolat,wet
  type(domain2D),       intent(in)    :: domain

  real, dimension(:),   allocatable          :: data
  real, dimension(:,:), allocatable          :: tmp
  real, dimension(:,:), allocatable          :: tmp_local
  real, dimension(:,:), allocatable          :: tmp1_local
  real, dimension(:,:), allocatable          :: tmp2_local
  real, dimension(:,:), allocatable          :: tmpx
  real, dimension(:,:), allocatable          :: tmpy
  
  integer :: i, j, k, lon_Tripol, ni, nj
  integer :: ioff, joff, iend, jend
  character(len=128) :: grd_file ! contains wet field
  character(len=128) :: topog_file
  
  ! set the grid points coordinates (degrees) and grid spacing (degrees)


  grd_file='INPUT/grid_spec.nc'
  topog_file='INPUT/topog.nc'
  
  ni = grid_ni;nj = grid_nj

  if (size(geolon,1).ne.ni+1) call mpp_error(FATAL,'array size mismatch in call to set_ocean_hgrid_arrays')
  if (size(geolon,2).ne.nj+1) call mpp_error(FATAL,'array size mismatch in call to set_ocean_hgrid_arrays')
  if (size(geolat,1).ne.ni+1) call mpp_error(FATAL,'array size mismatch in call to set_ocean_hgrid_arrays')
  if (size(geolat,2).ne.nj+1) call mpp_error(FATAL,'array size mismatch in call to set_ocean_hgrid_arrays')
 
  allocate (x_c(ni+1,nj+1))
  allocate (y_c(ni+1,nj+1))
  allocate (grid_x_c(ni))
  allocate (grid_y_c(nj))
  allocate (grid_wet(size(wet,1),size(wet,2)))

  !--- initialize grid data
  
  x_c=0.0;    y_c=0.0;     grid_x_c=0.0; grid_y_c=0.0 ; grid_wet=0.0

  select case( grid_version )
  case( VERSION_0 )
     call read_data(ocean_hgrid, "gridlon_c",      grid_x_c, no_domain = .true.)
     call read_data(ocean_hgrid, "gridlat_c",      grid_y_c, no_domain = .true.)
  case( VERSION_1 )
     call read_data(ocean_hgrid, "grid_x_C",      grid_x_c, no_domain = .true.)
     call read_data(ocean_hgrid, "grid_y_C",      grid_y_c, no_domain = .true.)
  case( VERSION_2 )
     allocate(tmpx(2*ni+1, 2*nj+1), tmpy(2*ni+1, 2*nj+1) )
     call read_data(ocean_hgrid, "x", tmpx, no_domain=.TRUE.)
     call read_data(ocean_hgrid, "y", tmpy, no_domain=.TRUE.)
     do i = 1, ni
        grid_x_c(i) = tmpx(2*i,  2) 
     enddo

     lon_Tripol = ni/4 ! 90 for 1 degree grid
     do j = 1, nj  
        grid_y_c(j) = tmpy(2*lon_Tripol+1 , 2*j)
     enddo
  end select

  allocate(tmp(ni+1,nj+1))
  !--- xc
  select case(grid_version)
  case(VERSION_0)
     call read_data(ocean_hgrid, 'geolon_c', tmp(1:ni,1:nj), no_domain = .true.)
  case(VERSION_1)
     call read_data(ocean_hgrid, 'x_C', tmp(1:ni,1:nj), no_domain = .true.)   
  case(VERSION_2)
     tmp(1:ni+1,1:nj+1) = tmpx(1:2*ni+1:2,1:2*nj+1:2)
  end select
  x_c = tmp

  !--- yc
  select case(grid_version)
  case(VERSION_0)
     call read_data(ocean_hgrid, 'geolat_c', tmp(1:ni,1:nj), no_domain = .true.)
  case(VERSION_1)
     call read_data(ocean_hgrid, 'y_C', tmp(1:ni,1:nj), no_domain = .true.)   
  case(VERSION_2)
     tmp(1:ni+1,1:nj+1) = tmpy(1:2*ni+1:2,1:2*nj+1:2)
  end select
  y_c = tmp

!--- wet / depth
   select case(grid_version)
     case(VERSION_0)
        call read_data(ocean_hgrid, 'wet', grid_wet, domain)
     case (VERSION_1)
        call read_data(ocean_hgrid, 'wet', grid_wet, domain)
     case (VERSION_2)
        if (field_exist(topog_file,'depth')) then
           call read_data(topog_file,'depth',grid_wet, domain)
           where ( grid_wet > 0.0 )
              grid_wet = 1.0
           elsewhere 
              grid_wet = 0.0
           end where
        else
           call mpp_error(FATAL,'depth field does not exist in topog_file')
        endif
   end select
    
   geolon=x_c; geolat=y_c; wet=grid_wet

end subroutine set_ocean_hgrid_arrays

end module ocean_grids_mod
