! TK: modified to run mixed layer experiments with either qflux restoring
!       and/or qflux adjustment based on a previous qflux restoring run.
!
!     Added do_restore_sst, do_qflux_adj namelist parameters.  
!
!     If do_restore_sst is true, read in the observed sst 
!        appropriate to current point in seasonal cycle. 
!        Compute the heat flux restoring term for the time period.
!        Add this term to the atm-ocean heat flux.
!
!     If do_qflux_adj is true, read in the qflux adjustment 
!        appropriate to current point in seasonal cycle. 
!        Add this term to the atm-ocean heat flux.
!        Updated name of qflux_adjustment input file 11/01/01.
!
!       Note: Avoid resetting the SST/seaice with climatology, since these
!           are now computed.  (i.e., spec_ice should be false)
!
!     TK  11/15/01
!     TK   3/07/02   Adapted to galway code, leaving out warm_restore code.
!     TK   7/31/02   Adapted to havana code (minor updates)
!     TK  11/25/02   Modify get_heat_flux_adjust routine to extract 
!                    data onto local grid rather than interpolate.  
!                    (Also requires some mods to update_ocean_model.)
!     TK   3/26/03   Modified branch code (qflux_tk ocean_mixed_layer)
!                    to be compatible with the default ocean_model.f90
!                    pre-Inchon "latest" version of ocean_model.f90 
!                    which uses data override.  That version was checked out
!                    using:  cvs co -r latest fms_spectral_mixed_layer.
!
!     TK   4/3/03    Added get_sea_surface routine directly to this
!                    module for use with the ocean_model.  Then
!                    data override will input data on the ocean
!                    grid with no conflict...
!
!     TK   4/14/03   Modified get_heat_flux_adj subroutine to work
!                    with data override.
!
!     TK   4/15/03   Updated model to be compatible with default
!                    inchon code.  
!                    Used inchon ocean_mixed_layer/ocean_model.f90
!                    as basis for this...
! 
!     TK   7/24/03   Code clean-up and documentation.  Subroutine
!                    get_heat_flux_adj was made not public.
!                    Added xml comments.  Removed extra mlcp definition
!                    statement.  Changed default mixed layer depth
!                    to 50 meters.  None of these should change answers
!                    and did not in a 1-year test run of the
!                    AM2p12.7_Res run (checksums).
!
!     TK   9/10/03   Added boolean module_is_initialized and 
!                    associated tests to bring up to FMS specs.
!                    Should not change answers.   
!
!     TK   1/09/04   moved get_sea_surface call inside if statement
!                    Should not change answers.
!
!     TK   3/10/04   added 5 more pointers related to shortwave
!                    radiation fluxes to the ice_ocean_boundary_type.
!                    This change originates from Rick Hemler's 
!                    updates to radiation code.
!
!     TK   4/27/04   added/modified code so that observed sst's
!                    can be saved out (i.e., the restoring target).
!
!     TK   8/31/04   modified subroutine ocean_diagnostics to produce
!                    netcdf output that is modulo in longitude by having
!                    edges information (xb,yb) for the x,y axes (to be
!                    consistent with ice model, and needed to work
!                    correctly in ferret, etc.)
!
!     LAT  4/4/07    split the sw_flux into 4 components to comply with 
!                    the Nalanda flux exchange; sw_flux_nir_dir,
!                    sw_flux_nir_dif, sw_flux_vis_dir, sw_flux_vis_dif. Removed
!                    sw_flux and sw_flux_vis.
!
!     LAT  5/12/08   modified the top-level ocean interfaces to Bob Hallberg's
!                    proposed quasi-object-oriented approach to facilitate such
!                    things as data assimilation, running with ensembles and 
!                    nesting. Time is passed around explicitly. ocean_data_type
!                    is renamed to ocean_public_type, with time, mask, and all 
!                    variables written to restart files moved to ocean_state_type,
!                    pe renamed to is_ocean_pe, and instance_name, BGRID_NE and
!                    avg_kount added. ocean_state_type is a new, private pointer containing
!                    time and mask. The collective 'data' field in ice_ocean_boundary_type
!                    was not used and is removed. ocean_model_init,
!                    update_ocean_model and ocean_model_end contain the
!                    additional ocean_state arugument. time_step was removed from 
!                    ocean_public_type and is an argument 'Ocean_coupling_time_step' in
!                    update_ocean_model, so it does not need to be stored. Ocean_seg_start,
!                    Ocean_seg_end and num_ocean_calls were unused place holders and removed
!                    as arguments from the update_ocean_model interface. The
!                    argument Ice_boundary was renamed to Ice_ocean_boundary. A new subroutine,
!                    ocean_model_save_restart, has been added for incremental
!                    checkpointing. These changes bitwise reproduce answers.
!
!     sdu  11/30/08  modified code to allow use of the mosiac gridspec type.
!
!     LAT  3/6/09    Added option to restore the SST to historical (time varying) data, rather
!                    than climatological data. At this time, the historical data is 
!                    HadISST1_SST_Dec1869-2004_1x1.nc. This is controlled via a namelist option,
!                    do_restore_histsst. For backward compatibility, do_restore_sst uses the
!                    climatological SST data. A fatal error will result if both of these
!                    namelists are set to true. The new diagnostic, sst_histobs, was created. 


module ocean_model_mod

  use    constants_mod, only:  Tfreeze, hlv, hlf, pi, rho0r   ! MS mod

  use time_manager_mod, only: time_type, operator(+), operator(>), get_date, &
                              get_time, set_time

  use mpp_mod,         only: mpp_pe, mpp_npes

  use mpp_domains_mod, only: mpp_domains_init, mpp_define_domains, domain2D, &
                             mpp_define_layout, mpp_get_compute_domain,      &
                             mpp_get_layout, FOLD_NORTH_EDGE, CYCLIC_GLOBAL_DOMAIN

  use mpp_parameter_mod, only: AGRID, BGRID_NE, CGRID_NE, BGRID_SW, CGRID_SW

  use fms_mod, only: close_file, check_nml_error, error_mesg, NOTE, &
                     WARNING, FATAL, write_version_number, stdlog, open_restart_file, &
                     file_exist, read_data, write_data, set_domain, mpp_pe,       &
                     mpp_root_pe, open_restart_file

  use diag_manager_mod, only: diag_axis_init, register_diag_field, send_data
#ifdef INTERNAL_FILE_NML
  USE mpp_mod, ONLY: input_nml_file
#else
  USE fms_mod, ONLY: open_namelist_file
#endif

  use ocean_grids_mod, only: ocean_grids_init, set_ocean_grid_size, set_ocean_hgrid_arrays

!   TK Mod: do not use ice_spec version of get_sea_surface in the
!      ocean model.  It has its own version (for data override).
!!  use ice_spec_mod,     only: get_sea_surface


! TK added this use statement, needed in local get_sea_surface routine...
  use data_override_mod,only: data_override
  use coupler_types_mod,only: coupler_2d_bc_type
! <CONTACT EMAIL="Tom.Knutson@noaa.gov">
!   Tom Knutson
! </CONTACT>

! <OVERVIEW>
!   This is a slab-type ocean model.  It is the primary means
!      of determining the equilibrium climate sensitivity 
!      of an atmosphere/land/sea-ice model.
! </OVERVIEW>

! <DESCRIPTION>
!   This is a slab-type ocean model, sometimes referred to as
!      as a qflux-adjusted mixed-layer model.  It is a very
!      simple ocean model consisting of a single motionless 
!      layer without advection or diffusion, but possibly 
!      including sea ice.  This type of model is the 
!      primary tool used to determine the equilibrium climate 
!      sensitivity of an atmosphere/land/sea-ice model.
!
!   The model is designed to run in 3 stages:  restoring
!      mode, qflux-adjusted control run mode, and perturbation
!      experiment mode.  The model can also be run in 
!      "free mode" with no restoring or qflux adjustments,
!      although this is normally not done.  It can also
!      be run in various other combinations such as qflux-adjusted
!      with further restoring although this is also not 
!      normally done.
!
!   Note that the sea-ice model is an integral part of the slab-ocean 
!      model although that code exists elsewhere (in the SIS part of 
!      FMS).  Care should be taken that the ice restoring fluxes 
!      archived by the sea-ice model during the restoring phase 
!      are added to the sst restoring fluxes archived by this
!      model during the restoring phase in order to determine 
!      the total qflux adjustment field.  Similarly any ice-lid fluxes
!      archived in the control experiment by the ice model
!      should be added to the control run qflux adjustment field
!      to obtain the final qflux adjustment field used 
!      for a perturbation experiment (without a lid).  More 
!      discussion of these topics is provided in the postscript 
!      documentation for the model (in preparation).
!
!
!   <LINK SRC="http://www.gfdl.noaa.gov/~tk/Slab_Model_Documentation.htm"> Click here for additional documentation </LINK>
! </DESCRIPTION>



  implicit none
  include 'netcdf.inc'
  private

!-----------------------------------------------------------------------
!----------------- public interfaces -----------------------------------

  public :: ocean_model_init, ocean_model_end, update_ocean_model,    &
            read_ice_ocean_boundary, write_ice_ocean_boundary, ocean_stock_pe,&
            init_default_ice_ocean_boundary, ocean_model_init_sfc,    &
            ocean_model_flux_init, ocean_model_save_restart, ocean_model_restart, &
            ocean_public_type_chksum, ice_ocn_bnd_type_chksum
!            get_heat_flux_adj        ! TK Mod (make not public) 7/24/03

public    ocean_model_data_get
interface ocean_model_data_get
   module procedure ocean_model_data1D_get 
   module procedure ocean_model_data2D_get 
end interface

!-----------------------------------------------------------------------
character(len=128), parameter :: version = '$Id$'
character(len=128), parameter :: tagname = '$Name$'

! <PUBLICTYPE>
  type, public :: ocean_public_type
     type(domain2D) :: Domain

     real, pointer, dimension(:,:) :: &
          t_surf =>NULL(), &      !  mixed layer temperature
          s_surf =>NULL(), &      !  mixed layer salinity (a constant at present in this model)
          sea_lev =>NULL(), &     !  sea level (Not used in Mixed layer model)
          frazil =>NULL(), &      !  energy flux of frazil formation.
                                  !     unit:  when divided by dt_ocean, unit is W/m2
          u_surf =>NULL(), &      !  zonal ocean current (= 0 at present in this model)
          area   =>NULL() , &
          v_surf =>NULL()         !  meridional ocean current (= 0 at present in this model)
     logical, pointer, dimension(:,:) :: maskmap =>NULL()! A pointer to an array indicating which
                                                         ! logical processors are actually used for
                                                         ! the ocean code. The other logical
                                                         ! processors would be all land points and
                                                         ! are not assigned to actual processors.
                                                         ! This need not be assigned if all logical
                                                         ! processors are used. This variable is dummy and need 
                                                         ! not to be set, but it is needed to pass compilation.
     integer, pointer, dimension(:) :: pelist =>NULL()
     logical :: is_ocean_pe           ! .true. on processors that run the ocean model.
     character(len=32) :: instance_name = ""   ! A name that can be used to identify
                                               ! this instance of an ocean model, for example
                                               ! in ensembles when writing messages.
     integer :: stagger = BGRID_NE    ! The staggering relative to the tracer
                                      ! points of the two velocity components.
                                      ! Valid entries include AGRID, BGRID_NE,
                                      ! CGRID_NE, BGRID_SW, and CGRID_SW, 
                                      ! corresponding to the community-standard
                                      ! Arakawa notation. (These are named 
                                      ! integers taken from mpp_parameter_mod.)
                                      ! Following MOM, this is BGRID_NE by
                                      ! default.
     integer, dimension(3)            :: axes = 0 ! Axis numbers that are
                                                  ! available for I/O using 
                                                  ! this surface data.  
     type(coupler_2d_bc_type)         :: fields   ! array of fields used for additional tracers
     integer                          :: avg_kount  ! Used for accumulating
                                                    ! averages of this type 
  end type ocean_public_type
! </PUBLICTYPE>

  type, public :: ocean_state_type ; private
     logical :: is_ocean_pe = .false.    ! .true. if this is an ocean PE.
     type(time_type) :: ocean_time       ! The ocean's internal clock - it is
                                         ! here for error checking against the
                                         ! input variables, and because defined
                                         ! types cannot be empty.
     logical, pointer, dimension(:,:) :: mask =>NULL() ! mask is true for ocean points, false for land points
  end type ocean_state_type


!Balaji

! <PUBLICTYPE>
  type, public :: ice_ocean_boundary_type      
     real, dimension(:,:), pointer :: u_flux =>NULL(), &    ! zonal wind stress (Pa)
                                      v_flux =>NULL(), &    ! meridional wind stress (Pa)
                                      t_flux =>NULL(), &    ! sensible heat flux (w/m2)
                                      q_flux =>NULL(), &    ! specific humidity flux (kg/m2/s)
                                      salt_flux =>NULL(), & ! salinity flux (kd/m2/s) (Not used in mixed layer model) 
                                      lw_flux =>NULL(), &   ! net (down-up) longwave flux (W/m2)
                                      sw_flux_nir_dir => NULL(),  & ! direct NIR sw radiation (W/m2) 
                                      sw_flux_nir_dif => NULL(),  & ! diffuse NIR sw radiation (W/m2)
                                      sw_flux_vis_dir => NULL(),  & ! direct visible sw radiation (W/m2)
                                      sw_flux_vis_dif => NULL(),  & ! diffuse visible sw radiation (W/m2)
                                      lprec =>NULL(), &     ! mass flux of liquid precipitation (Kg/m2/s)
                                      fprec =>NULL()        ! mass flux of frozen precipitation (Kg/m2/s)
     real, dimension(:,:), pointer :: runoff =>NULL(), &    ! mass flux of liquid runoff (Kg/m2/s)
                                      calving =>NULL()      ! mass flux of frozen runoff (Kg/m2/s)
     real, pointer, dimension(:,:) :: runoff_hflx     =>NULL() ! heat flux of liquid runoff (kg/m2/s) 
     real, pointer, dimension(:,:) :: calving_hflx    =>NULL() ! heat flux of frozen runoff (kg/m2/s) 
     real, dimension(:,:), pointer :: p =>NULL()   ! pressure on the surface of the ocean (Pa) 
                                                   ! (Not used in mixed layer model)
     real, pointer, dimension(:,:) :: mi              =>NULL() ! mass of overlying sea ice 
     integer                       :: wind_stagger = -999      ! member needed by coupler for SIS2. Not used.
     integer                       :: xtype   ! REGRID, REDIST or DIRECT
     type(coupler_2d_bc_type)      :: fluxes  ! array of fields used for additional tracers
     real, pointer, dimension(:,:) :: ustar_berg !> This is here for compatibility and is unused
     real, pointer, dimension(:,:) :: area_berg  !> This is here for compatibility and is unused
     real, pointer, dimension(:,:) :: mass_berg  !> This is here for compatibility and is unused
  end type ice_ocean_boundary_type
! </PUBLICTYPE>


  real, allocatable, dimension (:,:) :: geo_lonv, geo_latv, geo_lon, geo_lat
  real, allocatable, dimension (:,:  ) :: dum1
  real, allocatable, dimension (:,:,:) :: dum2
  integer :: im, jm
  real, allocatable, dimension(:) :: xb1d, yb1d ! 1d global grid for diag_manager

  integer :: id_sst, id_hflx, id_swflx_nir_dir, id_swflx_nir_dif, id_swflx_vis_dir, id_swflx_vis_dif, &
             id_lwflx, id_shflx, id_lhflx, id_snwflx, id_frazil       
  integer :: id_qflx_adj, id_qflux_restore_sst, id_net_hflx     ! TK Mod
  integer :: id_sst_obs     ! TK mod
  integer :: id_sst_histobs ! LAT mod: 3/6/09 - historical SST diagnostic
  integer :: id_pme, id_river     ! MS mod

  logical :: sent
  real    :: mlcp
  logical :: module_is_initialized = .FALSE.
  logical :: stock_warning_issued  = .FALSE.
!
! namelist
!
  real    :: mixed_layer_depth = 50.
  real    :: mixed_layer_salin = 33.333
  integer :: layout(2)=(/0,0/)
  logical :: do_qflux_adj     = .false.       ! TK mod
  logical :: do_restore_sst   = .false.       ! TK mod
  logical :: do_restore_histsst = .false.     ! LAT mod: 3/6/09 historical SST
  real    :: sst_restore_timescale   =  5.0   ! TK mod; restoring time scale for SST in days
  real    :: uniform_init_t_surf = -273.0
! <NAMELIST NAME="ocean_model_nml">
!   <DATA NAME="mixed_layer_depth" UNITS="meters"
!         TYPE="real" DEFAULT="50.0">
!     Depth of the mixed layer
!   </DATA>
!   <DATA NAME="mixed_layer_salin" UNITS="parts per thousand"
!         TYPE="real" DEFAULT="33.333">
!     Salinitiy of the mixed layer
!   </DATA>
!   <DATA NAME="layout" UNITS=" "
!         TYPE="integer" DEFAULT="(/0,0/)">
!     Multiple processor layout variable
!   </DATA>
!   <DATA NAME="do_qflux_adj" UNITS=" "
!         TYPE="logical" DEFAULT=".false.">
!     Add qflux adjustment to mixed layer tendency equation?
!   </DATA>
!   <DATA NAME="do_restore_sst" UNITS=" "
!         TYPE="logical" DEFAULT=".false.">
!     Restore mixed layer temperature toward observed SST?
!   </DATA>
!   <DATA NAME="do_restore_histsst" UNITS=" "
!         TYPE="logical" DEFAULT=".false.">
!     Restore mixed layer temperature toward observed historical SST? In this
!     case, the HadISST1_SST_Dec1869-2004_1x1.nc data.
!   </DATA>
!   <DATA NAME="sst_restore_timescale" UNITS="days"
!         TYPE="real" DEFAULT="5.0.">
!     Time scale for mixed layer temperature restoring (if used)
!   </DATA>
! </NAMELIST>

  namelist /ocean_model_nml/ mixed_layer_depth, mixed_layer_salin, layout, &
                             do_qflux_adj, do_restore_sst, do_restore_histsst,&
                             sst_restore_timescale, uniform_init_t_surf 
                             ! LAT Mod: 3/6/09 - added do_restore_histsst

!-----------------------------------------------------------------------
!-------- mpp domain2d types ---------

contains

! <SUBROUTINE NAME="ocean_model_init">

!   <OVERVIEW>
!     This initializes the ocean model
!   </OVERVIEW>

!   <DESCRIPTION>
!     The ocean model is initialized by reading namelist, 
!      reading in restart data, setting up variables and domains,
!      etc.
!   </DESCRIPTION>


! <PUBLICROUTINE>
  subroutine ocean_model_init (Ocean_sfc, Ocean_state, Time_init, Time_in)

    type (ocean_public_type), target, &
         intent(inout) :: Ocean_sfc   ! Ocean definition variable
    type (ocean_state_type), pointer  &
                       :: Ocean_state ! Ocean internal state
    type (time_type), intent(in) :: &
         Time_init, &               ! Initial time (currently not used by routine)
         Time_in                    ! Current time - the time at which to initialize
                                    ! the ocean model.
! </PUBLICROUTINE>

    real, allocatable, dimension(:,:) :: rmask
    integer :: unit, ierr, io
    integer :: i, j, is, ie, js, je, domain_flags
    integer :: rcode, ncid, varid, dims(4), start(4), nread(4)
    character(len=80) :: domainname

    module_is_initialized = .TRUE.
    
    allocate (Ocean_state)
    Ocean_state%is_ocean_pe = Ocean_sfc%is_ocean_pe
    if (.not.ocean_state%is_ocean_pe) return
    Ocean_state%ocean_time  = Time_in
!
! read namelist
!
#ifdef INTERNAL_FILE_NML
    READ (input_nml_file, NML=ocean_model_nml, IOSTAT=io)
    IF ( io > 0 ) CALL error_mesg('ocean_model_mod',&
         & 'Error reading internal namelist ocean_model_nml.', FATAL)
#else
    unit = open_namelist_file()
    ierr=1
    do while (ierr /= 0)
      read(unit, nml=ocean_model_nml, iostat=io, end=10)
      ierr = check_nml_error (io, 'ocean_model_nml')
    enddo
10  call close_file (unit)
#endif
!
! write version number and namelist
!
    call write_version_number(version, tagname)
    unit = stdlog()
    if ( mpp_pe() == mpp_root_pe() ) write (unit, nml=ocean_model_nml)

! Initilize the grid
    call ocean_grids_init()
    call set_ocean_grid_size(im, jm, 'INPUT/grid_spec.nc')

!
! set up domains
!
    call mpp_domains_init
    if( layout(1).EQ.0 .AND. layout(2).EQ.0 )call mpp_define_layout( (/1,im,1,jm/), mpp_npes(), layout )
    if( layout(1).NE.0 .AND. layout(2).EQ.0 )layout(2) = mpp_npes()/layout(1)
    if( layout(1).EQ.0 .AND. layout(2).NE.0 )layout(1) = mpp_npes()/layout(2)
    domainname = 'ML Ocean'
    call mpp_define_domains( (/1,im,1,jm/), layout, Ocean_sfc%Domain, xflags=CYCLIC_GLOBAL_DOMAIN, &
         name=domainname, yflags=FOLD_NORTH_EDGE )
!         yflags=CYCLIC_GLOBAL_DOMAIN )
    call mpp_get_compute_domain( Ocean_sfc%Domain, is, ie, js, je )

    allocate ( rmask (is:ie, js:je), Ocean_state%mask (is:ie, js:je) )
    allocate ( geo_lonv (1:im+1,1:jm+1), geo_latv (1:im+1,1:jm+1) )

    call set_ocean_hgrid_arrays(geo_lonv, geo_latv, rmask, Ocean_sfc%Domain)

    Ocean_state%mask = .false.; where (rmask > 0.5) Ocean_state%mask = .true.
    deallocate ( rmask )

    allocate ( geo_lon  (is:ie  , js:je  ), geo_lat  (is:ie  , js:je  ) )

    geo_lon = (geo_lonv(is:ie,js:je)    +geo_lonv(is+1:ie+1,js:je) &             
         +geo_lonv(is:ie,js+1:je+1)+geo_lonv(is+1:ie+1,js+1:je+1))/4
    geo_lat = (geo_latv(is:ie,js:je)    +geo_latv(is+1:ie+1,js:je) &
         +geo_latv(is:ie,js+1:je+1)+geo_latv(is+1:ie+1,js+1:je+1))/4

    allocate ( xb1d (im+1) )
    allocate ( yb1d (jm+1) )
    xb1d = sum(geo_lonv,2)/(jm+1);
    yb1d = sum(geo_latv,1)/(im+1);

    deallocate (geo_lonv)
    deallocate (geo_latv)

    mlcp = mixed_layer_depth*4e6

    allocate( Ocean_sfc%t_surf (is:ie,js:je), &
              Ocean_sfc%s_surf (is:ie,js:je), &
              Ocean_sfc%u_surf (is:ie,js:je), &
              Ocean_sfc%v_surf (is:ie,js:je), &
              Ocean_sfc%sea_lev(is:ie,js:je), &
              Ocean_sfc%frazil (is:ie,js:je)  )    

    Ocean_sfc%sea_lev = 0.
    Ocean_sfc%s_surf  = mixed_layer_salin;

    allocate ( dum1(is:ie,js:je), dum2(is:ie,js:je,2) )

!  The call to get_sea_surface is commented out because get_sea_surface
!  is also called in the sis ice model and it cannot put the fields on
!  the correct domain decomposition for both the ocean and ice models
!  at the same time. This is because get_sea_surface calls data_override,
!  which puts each override field on only one model domain, which is
!  specified in the data_table.
!  
!  A getaround is to initialize ssts to a simple function of latitude.

!   call get_sea_surface(Ocean_state%ocean_time, Ocean_sfc%t_surf, dum2, dum1)

    if ( file_exist('INPUT/ocean_model.res')) then

!         ---- set domain for global i/o ----
        call set_domain ( Ocean_sfc%Domain )

        unit = open_restart_file ('INPUT/ocean_model.res','read')

        call read_data ( unit, Ocean_sfc%t_surf )
        call read_data ( unit, Ocean_sfc%u_surf )
        call read_data ( unit, Ocean_sfc%v_surf )
        call read_data ( unit, Ocean_sfc%frazil )
        call close_file (unit)
    else
!     Cold start ssts are 0C poleward of 65 deg lat and 28C at equator
        Ocean_sfc%t_surf = 34*cos(pi*geo_lat/180.)**2 - 6. + Tfreeze
        if(uniform_init_t_surf > -100.0 ) Ocean_sfc%t_surf = Tfreeze + uniform_init_t_surf  
        where (Ocean_sfc%t_surf < Tfreeze) Ocean_sfc%t_surf = Tfreeze
        Ocean_sfc%u_surf  = 0.
        Ocean_sfc%v_surf  = 0.
        Ocean_sfc%frazil  = 0.
    endif

    call ocean_diagnostics (Ocean_sfc, Ocean_state%ocean_time)

  end subroutine ocean_model_init
! </SUBROUTINE>


! <SUBROUTINE NAME="update_ocean_model">

!   <OVERVIEW>
!     The main time-stepping routine for the mixed layer ocean model
!   </OVERVIEW>

!   <DESCRIPTION>
!     The mixed layer ocean model is integrated forward for
!      one ocean time step according to the following key set of 
!      equations:
!
!     hflx = Ice_ocean_boundary%sw_flux + Ice_ocean_boundary%lw_flux - &
!             (Ice_ocean_boundary%fprec+Ice_ocean_boundary%calving)*hlf - &
!             Ice_ocean_boundary%t_flux - Ice_ocean_boundary%q_flux*hlv
!
!     net_hflx(:,:) = hflx(:,:) + qflux_adj(:,:) + qflux_restore_sst(:,:)
!
!     Ocean_sfc%t_surf = Ocean_sfc%t_surf + dt_ocean * net_hflx/mlcp
!
!   </DESCRIPTION>


! <PUBLICROUTINE>
  subroutine update_ocean_model( Ice_ocean_boundary, Ocean_state, Ocean_sfc, &
       time_start_update, Ocean_coupling_time_step )

    type (ice_ocean_boundary_type), intent(inout) :: &
               Ice_ocean_boundary  ! Variable with fields for communication 
                                   ! from ice to the ocean
    type (ocean_state_type), pointer              :: &
               Ocean_state         ! Internal ocean state
    type (ocean_public_type), intent(inout) ::  &
               Ocean_sfc           ! Ocean definition variable 
    type (time_type), intent(in)          ::  &
               time_start_update,             &  ! time at beginning of update step
               Ocean_coupling_time_step          ! amount of time over which to advance ocean
! </PUBLICROUTINE>

    integer :: is, ie, js, je

    real, dimension(size(Ocean_state%mask,1),size(Ocean_state%mask,2)) :: sst, hflx
    real, dimension(size(Ocean_state%mask,1),size(Ocean_state%mask,2)) :: qflux_restore_sst, & ! TK mod  
                                                              qflux_adj, net_hflx  ! TK Mod
    real, dimension(size(Ocean_state%mask,1),size(Ocean_state%mask,2)) :: pme, river           ! MS mod
    real    :: dt_ocean
    integer :: dy, sc
    
      Ocean_state%ocean_time = Ocean_state%ocean_time + Ocean_coupling_time_step
      call get_time(Ocean_coupling_time_step,sc,dy)
      dt_ocean = 86400*dy + sc

!      if (time_start_update /= Ocean_state%ocean_time) then
!         call error_mesg('update_ocean_model', 'internal clock does not agree with time_start_update argument.', WARNING)
!      endif

    !Check if ocean_model has been initialized:
    if ( .NOT. module_is_initialized ) &
      !<ERROR MSG="Slab ocean_model not correctly initialized." STATUS="FATAL">
      ! The slab or mixed-layer type ocean model was not correctly initialized.
      ! </ERROR>
        call error_mesg( 'update_ocean_model', 'Slab ocean_model not correctly initialized.', FATAL )

      !check if required boundary fields have been initialized
      if( .NOT.ASSOCIATED(Ice_ocean_boundary%sw_flux_nir_dir) .OR. &
          .NOT.ASSOCIATED(Ice_ocean_boundary%sw_flux_nir_dif) .OR. &
          .NOT.ASSOCIATED(Ice_ocean_boundary%sw_flux_vis_dir) .OR. &
          .NOT.ASSOCIATED(Ice_ocean_boundary%sw_flux_vis_dif) .OR. &
          .NOT.ASSOCIATED(Ice_ocean_boundary%lw_flux) .OR. &
          .NOT.ASSOCIATED(Ice_ocean_boundary%fprec)   .OR. &
          .NOT.ASSOCIATED(Ice_ocean_boundary%calving) .OR. &
          .NOT.ASSOCIATED(Ice_ocean_boundary%t_flux)  .OR. &
          .NOT.ASSOCIATED(Ice_ocean_boundary%q_flux) ) &
      !<ERROR MSG="Ice_ocean_boundary not correctly initialized." STATUS="FATAL">
      ! The Ice_ocean_boundary variable was not correctly initialized.
      !</ERROR>
         call error_mesg( 'Update_ocean_model', 'Ice_ocean_boundary not correctly initialized.', FATAL )

      pme = 0.0        ! MS mod
      river = 0.0      ! MS mod
      hflx = 0.0
      net_hflx = 0.0   ! TK Mod
      ! mlcp = 4e6*mixed_layer_depth   TK Mod: removed extra statement.
      ! TK Mods to implement qflux adjusted or restoring mixed layer models
      ! Original Galway Code (next 5 lines) for reference with SW flux partition from Nalanda
      ! where (Ocean_state%mask)
      !   hflx = Ice_ocean_boundary%sw_flux_nir_dir + Ice_ocean_boundary%sw_flux_nir_dif + &
      !   Ice_ocean_boundary%sw_flux_vis_dir + Ice_ocean_boundary%sw_flux_vis_dif + &
      !   Ice_ocean_boundary%lw_flux - (Ice_ocean_boundary%fprec+Ice_ocean_boundary%calving)*hlf - &
      !   Ice_ocean_boundary%t_flux - Ice_ocean_boundary%q_flux*hlv 
      !   Ocean_sfc%t_surf = Ocean_sfc%t_surf + dt_ocean*hflx/mlcp
      ! end where


      qflux_restore_sst(:,:) = 0.0

! LAT mod: 3/6/09 - SST restoring must be either from the climatology or 
! historical datasets.
      if (do_restore_sst .and. do_restore_histsst) then
         call error_mesg('Update_ocean_model', &
                     'do_restore_sst and do_restore_histsst cannot both be set to true.', FATAL)
      endif
 
! LAT mod: 3/6/09 - added option for restoring to historical SSTs
      if (do_restore_sst .or. do_restore_histsst) then
     
! TK Mod (1/09/2004) moved this statement inside this if statement
!    since it is only needed if sst restoring is active.
!    (for spec_ice = true, the sst's are read in in ice_model.f90)  

        call get_sea_surface(Ocean_state%ocean_time, sst, dum2, dum1)

!  TK Mod: 4/24/04
        if (id_sst_obs    > 0) &
          sent = send_data(id_sst_obs, sst, Ocean_state%ocean_time, mask=Ocean_state%mask)
!  LAT Mod: 3/6/09 - added diagnostic id for historical restoring SSTs
        if (id_sst_histobs    > 0) &
          sent = send_data(id_sst_histobs, sst, Ocean_state%ocean_time, mask=Ocean_state%mask)
      
          ! Find the heat flux restoring term (qflux_restore_sst) due to sst
          ! restoring. The qflux_restore_sst values are saved for use in later
          ! flux adjusted runs. To find the sst restoring, use the observed sst
          ! for the date (stored in sst, which has a dual role in this routine)
    
          where (Ocean_state%mask)
            qflux_restore_sst = (sst - Ocean_sfc%t_surf) &
                                 * mlcp / (sst_restore_timescale * 86400)                  
          end where
      end if    

      if (do_qflux_adj) then
  
  !   Find the current climatological heat flux adjustment (qflux_adj).  
  !   The qflux_adj can be added to the heat flux to find the net
  !   heat flux anomaly into the mixed layer for the coupled mixed layer expt.
  
  !   Note that any ice model contributions need to be added in, which
  !   is handled in an off-line calculation to derive the final qflux
  !   adjustment.  These include ice restoring and ice lid contributions.

  !   Obtain indices for local array:
        call mpp_get_compute_domain( Ocean_sfc%Domain, is, ie, js, je )

        call get_heat_flux_adj (Ocean_state%ocean_time, qflux_adj)       

      else
        qflux_adj = 0.0
      end if    
    
      where (Ocean_state%mask)

! <PUBLICCOMMENT>
! The following are some key equations for this model:
     
       hflx = Ice_ocean_boundary%sw_flux_nir_dir + Ice_ocean_boundary%sw_flux_nir_dif + &
              Ice_ocean_boundary%sw_flux_vis_dir + Ice_ocean_boundary%sw_flux_vis_dif + &
              Ice_ocean_boundary%lw_flux - (Ice_ocean_boundary%fprec+Ice_ocean_boundary%calving)*hlf - &
              Ice_ocean_boundary%t_flux - Ice_ocean_boundary%q_flux*hlv

       net_hflx(:,:) = hflx(:,:) + qflux_adj(:,:) + qflux_restore_sst(:,:)

       Ocean_sfc%t_surf = Ocean_sfc%t_surf + dt_ocean*net_hflx/mlcp
     
       pme = (Ice_ocean_boundary%lprec + Ice_ocean_boundary%fprec - Ice_ocean_boundary%q_flux)*rho0r   ! MS mod
       river = (Ice_ocean_boundary%runoff + Ice_ocean_boundary%calving)*rho0r                    ! MS mod

! </PUBLICCOMMENT>     
      end where

!   End of TK Mods Section

      where (Ocean_state%mask .and. Ocean_sfc%t_surf < Tfreeze-0.054*mixed_layer_salin)
          Ocean_sfc%frazil = (Tfreeze-0.054*mixed_layer_salin-Ocean_sfc%t_surf)*mlcp
          Ocean_sfc%t_surf =  Tfreeze-0.054*mixed_layer_salin
      else where
          Ocean_sfc%frazil = 0.0
      end where

!-----------------------------------------------------------------------

      sst = 0.0; where (Ocean_state%mask) sst = Ocean_sfc%t_surf-Tfreeze

      if (id_sst    > 0) &
           sent = send_data(id_sst, sst, Ocean_state%ocean_time, mask=Ocean_state%mask)
      if (id_hflx   > 0) &
           sent = send_data(id_hflx, hflx, Ocean_state%ocean_time, mask=Ocean_state%mask)
      if (id_swflx_nir_dir  > 0) &
           sent = send_data(id_swflx_nir_dir, Ice_ocean_boundary%sw_flux_nir_dir, Ocean_state%ocean_time, mask=Ocean_state%mask)
      if (id_swflx_nir_dif  > 0) &
           sent = send_data(id_swflx_nir_dif, Ice_ocean_boundary%sw_flux_nir_dif, Ocean_state%ocean_time, mask=Ocean_state%mask)
      if (id_swflx_vis_dir  > 0) &
           sent = send_data(id_swflx_vis_dir, Ice_ocean_boundary%sw_flux_vis_dir, Ocean_state%ocean_time, mask=Ocean_state%mask)
      if (id_swflx_vis_dif  > 0) &
           sent = send_data(id_swflx_vis_dif, Ice_ocean_boundary%sw_flux_vis_dif, Ocean_state%ocean_time, mask=Ocean_state%mask)
      if (id_lwflx  > 0) &
           sent = send_data(id_lwflx, Ice_ocean_boundary%lw_flux, Ocean_state%ocean_time, mask=Ocean_state%mask)
      if (id_shflx  > 0) &
           sent = send_data(id_shflx, -Ice_ocean_boundary%t_flux, Ocean_state%ocean_time, mask=Ocean_state%mask)
      if (id_lhflx  > 0) &
           sent = send_data(id_lhflx, -Ice_ocean_boundary%q_flux*hlv, Ocean_state%ocean_time, mask=Ocean_state%mask)
      if (id_snwflx > 0) &
           sent = send_data(id_snwflx, -Ice_ocean_boundary%fprec*hlf, Ocean_state%ocean_time, mask=Ocean_state%mask)
      if (id_frazil > 0) &
           sent = send_data(id_frazil, Ocean_sfc%frazil/dt_ocean, Ocean_state%ocean_time, &
           mask=Ocean_state%mask)
! MS mods
      if (id_pme    > 0) &
           sent = send_data(id_pme, pme, Ocean_state%ocean_time, mask=Ocean_state%mask)
      if (id_river  > 0) &
           sent = send_data(id_river, river, Ocean_state%ocean_time, mask=Ocean_state%mask)

! TK mods
      if (do_qflux_adj) then
        if (id_qflx_adj   > 0) &
           sent = send_data(id_qflx_adj, qflux_adj, Ocean_state%ocean_time, mask=Ocean_state%mask)
      end if      
   
! LAT mod: 3/6/09 - added option to send data for restoring historical SSTs
      if (do_restore_sst .or. do_restore_histsst) then
        if (id_qflux_restore_sst   > 0) &
           sent = send_data(id_qflux_restore_sst, qflux_restore_sst, Ocean_state%ocean_time, mask=Ocean_state%mask)
      end if      
   
      if (do_qflux_adj .or. do_restore_sst .or. do_restore_histsst) then
        if (id_net_hflx   > 0) &
           sent = send_data(id_net_hflx, net_hflx, Ocean_state%ocean_time, mask=Ocean_state%mask)
      end if      

  end subroutine update_ocean_model
! </SUBROUTINE>

! <SUBROUTINE NAME="ocean_model_end">

!   <OVERVIEW>
!     End of integration routine for the mixed layer ocean model
!   </OVERVIEW>

!   <DESCRIPTION>
!     This routine writes and closes the restart file and is
!     called at the end of a model run on the system.

!   </DESCRIPTION>


! <PUBLICROUTINE>
  subroutine ocean_model_end (Ocean_sfc, Ocean_state, Time)

    type (ocean_public_type),           intent(inout) :: &
                     Ocean_sfc        ! Ocean definition variable
    type (ocean_state_type),            pointer       :: &
                     Ocean_state      ! Internal ocean state
    type (time_type),                   intent(in)    :: &
                     Time             ! Model time, used for writing restarts
! </PUBLICROUTINE>
    integer :: yr, mo, dy, hr, mn, sc
    integer :: unit

!    ---- set domain for global i/o ----
    call set_domain ( Ocean_sfc%Domain )

    unit = open_restart_file ('RESTART/ocean_model.res','write')

    call write_data ( unit, Ocean_sfc%t_surf )
    call write_data ( unit, Ocean_sfc%u_surf )
    call write_data ( unit, Ocean_sfc%v_surf )
    call write_data ( unit, Ocean_sfc%frazil )

    call close_file (unit)

    module_is_initialized = .FALSE.
    
  end subroutine ocean_model_end
! </SUBROUTINE>
 
!#######################################################################
! <SUBROUTINE NAME="ocean_model_restart">
!
! <DESCRIPTION>
! dummy interface.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 
! </DESCRIPTION>
!
  subroutine ocean_model_restart(Ocean_state, timestamp)
     type(ocean_state_type),    pointer     :: Ocean_state
     character(len=*), intent(in), optional :: timestamp

    call error_mesg('ocean_model_restart(ocean_model_mod)', &
                     'intermediate restart capability is not implemented for this model', FATAL)

  end subroutine ocean_model_restart
! </SUBROUTINE> NAME="ocean_model_restart"

! <SUBROUTINE NAME="ocean_model_save_restart">

!   <OVERVIEW>
!     Routine to save ocean restart files.
!   </OVERVIEW>

!   <DESCRIPTION>
!     This routine writes restart data to restart files.
!     It exists for incremental checkpointing, as well as
!     writing restart files at the end of the run from 
!     ocean_model_end.
!   </DESCRIPTION>


! <PUBLICROUTINE>
  subroutine ocean_model_save_restart(Ocean_state, Time, directory, filename_suffix)
    type(ocean_state_type),     pointer    :: &
              Ocean_state     ! Internal ocean state
    type(time_type),            intent(in) :: &
              Time            ! Model time at this call. Needed for mpp_write calls.
    character(len=*), optional, intent(in) :: &
              directory       ! Optional directory into which to write restart files
    character(len=*), optional, intent(in) :: &
              filename_suffix ! Optional suffix (e.g., a time-stamp) to append
                              ! to the restart file names
! </PUBLICROUTINE>
    return
  end subroutine ocean_model_save_restart
! </SUBROUTINE>

  
  subroutine ocean_diagnostics (Ocean_sfc, Time)
    type (ocean_public_type), intent(in) :: Ocean_sfc
    type (time_type), intent(in) :: Time

! TK mods on 8/30/04:  needed to produce modulo netcdf output for ocean model.
    integer, dimension(2) :: axt    
    integer id_xb, id_xt, id_yb, id_yt

!     id_xt = diag_axis_init('xto', (xb1d(1:im)+xb1d(2:im+1))/2, 'degrees_E', 'x', &
!          set_name='ocean', Domain2=Ocean_sfc%Domain )
!     id_yt = diag_axis_init('yto', (yb1d(1:jm)+yb1d(2:jm+1))/2, 'degrees_N', 'y', &
!          set_name='ocean', Domain2=Ocean_sfc%Domain )
     id_xb = diag_axis_init('xb', xb1d, 'degrees_E', 'X', 'longitude', &
                                    set_name='ocean',Domain2=Ocean_sfc%Domain )
     id_xt = diag_axis_init('xto', (xb1d(1:im)+xb1d(2:im+1))/2, 'degrees_E', 'x', &
          'longitude',set_name='ocean',edges=id_xb,Domain2=Ocean_sfc%Domain )
     id_yb = diag_axis_init('yb', yb1d, 'degrees_N', 'Y', 'latitude', &
                                   set_name='ocean',Domain2=Ocean_sfc%Domain )
     id_yt = diag_axis_init('yto', (yb1d(1:jm)+yb1d(2:jm+1))/2, 'degrees_N', 'y', &
          'latitude',set_name='ocean',edges=id_yb,Domain2=Ocean_sfc%Domain )

    axt  = (/ id_xt, id_yt       /);

!    Sample of old style:  (changed (/id_xt,id_yt/) to axt)
!    id_sst = register_diag_field('ocean_model','SST',(/id_xt,id_yt/), Time, &
!         'sea surface temp.',  'deg-C', missing_value=1.e10)

    id_sst = register_diag_field('ocean_model','SST', axt, Time, &
         'sea surface temp.',  'deg-C', missing_value=1.e10)

    id_hflx = register_diag_field('ocean_model','HFLX',axt,Time,&
         'surface heat flux',  'W/m^2', missing_value=1.e10)

    id_swflx_nir_dir = register_diag_field('ocean_model','SWFLX_nir_dir',axt, &
         Time, 'short wave NIR direct flux', 'W/m^2', missing_value=1.e10)

    id_swflx_nir_dif = register_diag_field('ocean_model','SWFLX_nir_dif',axt, &
         Time, 'short wave NIR diffuse flux', 'W/m^2', missing_value=1.e10)

    id_swflx_vis_dir = register_diag_field('ocean_model','SWFLX_vis_dir',axt, &
         Time, 'short wave visible direct flux', 'W/m^2', missing_value=1.e10)

    id_swflx_vis_dif = register_diag_field('ocean_model','SWFLX_vis_dif',axt, &
         Time, 'short wave visible diffuse flux', 'W/m^2', missing_value=1.e10)

    id_lwflx = register_diag_field('ocean_model','LWFLX',axt, &
         Time, 'long wave flux', 'W/m^2', &
         missing_value=1.e10)

    id_shflx = register_diag_field('ocean_model','SHFLX',axt, &
         Time, 'sensible heat flux', 'W/m^2', &
         missing_value=1.e10)

    id_lhflx = register_diag_field('ocean_model','LHFLX',axt, &
         Time, 'latent heat flux','W/m^2', &
         missing_value=1.e10)

    id_snwflx = register_diag_field('ocean_model','SNWFLX',axt, &
         Time,'latent heat of snow','W/m^2', &
         missing_value=1.e10)

    id_frazil = register_diag_field('ocean_model','FRAZIL'  ,axt, &
         Time, 'frazil', 'W/m^2', &
         missing_value=1.e10)

!   MS mods
    id_pme = register_diag_field('ocean_model','PME', axt, Time, &
         'prec-evap (liquid, frozen, evapor)',  'm/sec', missing_value=1.e10)

    id_river = register_diag_field('ocean_model','RIVER', axt, Time, &
         'river flux (liquid, frozen)',  'm/sec', missing_value=1.e10)

!   TK mods
!   LAT mod: 3/6/09 - added calls/diagnostics for restoring historical SST case
    if (do_restore_sst .or. do_restore_histsst)  then
         id_qflux_restore_sst = register_diag_field('ocean_model', &
         'QFLX_RESTORE_SST',axt, &
         Time, 'surface heat flux for SST restoring',  'W/m^2', &
         missing_value=1.e10)
         if (do_restore_histsst) then
            call error_mesg ('ocean_model_mod', 'restoring to historical SSTs', NOTE)
            id_sst_histobs = register_diag_field('ocean_model', &
            'SST_HISTOBS',axt, &
            Time, 'Observed historical SST for restoring',  'K', &
            missing_value=1.e10)  
         else
            call error_mesg ('ocean_model_mod', 'restoring to climatological SSTs', NOTE)
            id_sst_obs = register_diag_field('ocean_model', &
            'SST_OBS',axt, &
            Time, 'Observed SST for restoring',  'K', &
            missing_value=1.e10)
         endif
    end if      

    if (do_qflux_adj) then
         id_qflx_adj = register_diag_field('ocean_model', &
         'QFLUX_ADJ',axt, &
         Time, 'surface heat flux adjustment',  'W/m^2', &
         missing_value=1.e10)
    end if                          

! LAT mod: 3/6/09 - added restoring historical SST case
    if (do_qflux_adj .or. do_restore_sst .or. do_restore_histsst) then
         id_net_hflx = register_diag_field('ocean_model', &
         'NET_HFLX',axt, &
         Time, 'net heat flux into mixed layer after adjustments',  &
         'W/m^2', missing_value=1.e10)
    end if                          

  end subroutine ocean_diagnostics

! <SUBROUTINE NAME="read_ice_ocean_boundary">

!   <OVERVIEW>
! This ocean_model doesn't support concurrent runs,
! so this routine should never be called. This is a dummy routine.
!   </OVERVIEW>

!   <DESCRIPTION>
! This ocean_model doesn't support concurrent runs,
! so this routine should never be called. This is a dummy routine.

!   </DESCRIPTION>


! <PUBLICROUTINE>

  subroutine read_ice_ocean_boundary(file_name,iob,Ocean)

    character(LEN=*),             intent(IN)    :: file_name  
    type(ice_ocean_boundary_type),intent(INOUT) :: iob
    type(ocean_public_type),      intent(IN)    :: Ocean
! </PUBLICROUTINE>

!   <ERROR MSG="ocean_mixed_layer model cannot run concurrently." STATUS="FATAL">
!     The ocean_mixed_layer model is not coded to run concurrently at present.
!   </ERROR>
    call error_mesg('ocean_model_mod', 'ocean_mixed_layer model cannot run concurrently', &
                    FATAL )
    return

  end subroutine read_ice_ocean_boundary
! </SUBROUTINE>

! <SUBROUTINE NAME="read_ice_ocean_boundary">

!   <OVERVIEW>
! This ocean_model doesn't support concurrent runs,
! so this routine should never be called. This is a dummy routine.
!   </OVERVIEW>

!   <DESCRIPTION>
! This ocean_model doesn't support concurrent runs,
! so this routine should never be called. This is a dummy routine.

!   </DESCRIPTION>


! <PUBLICROUTINE>

  subroutine write_ice_ocean_boundary(file_name,iob,Ocean)

    character(LEN=*),             intent(IN) :: file_name  
    type(ice_ocean_boundary_type),intent(IN) :: iob
    type(ocean_public_type),      intent(IN) :: Ocean
! </PUBLICROUTINE>

!   <ERROR MSG="ocean_mixed_layer model cannot run concurrently." STATUS="FATAL">
!     The ocean_mixed_layer model is not coded to run concurrently at present.
!   </ERROR>
    call error_mesg('ocean_model_mod', 'ocean_mixed_layer model cannot run concurrently', &
                    FATAL )
    return

  end subroutine write_ice_ocean_boundary
! </SUBROUTINE>

! <SUBROUTINE NAME="init_default_ice_ocean_boundary">

!   <OVERVIEW>
! dummy routine
!   </OVERVIEW>

!   <DESCRIPTION>
! dummy routine
!   </DESCRIPTION>


! <PUBLICROUTINE>

! dummy routine

  subroutine init_default_ice_ocean_boundary(iob)

    type(ice_ocean_boundary_type),intent(INOUT) :: iob
! </PUBLICROUTINE>

    return

  end subroutine init_default_ice_ocean_boundary
! </SUBROUTINE>

! dummy interface for ESM coupler
subroutine ocean_model_init_sfc(Ocean_state,Ocean)

    type(ocean_state_type),  pointer             :: Ocean_state
    type(ocean_public_type), intent(inout)       :: Ocean
    
return
end subroutine ocean_model_init_sfc

subroutine ocean_model_flux_init(Ocean_state)
  type(ocean_state_type), pointer :: Ocean_state  ! Internal ocean state

  if (.not.associated(Ocean_state)) then
    allocate(Ocean_state)
    Ocean_state%is_ocean_pe = .false.
  endif
return
end subroutine ocean_model_flux_init


!#######################################################################



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! get_sea_surface - get SST, ice concentration and thickness from data         !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine get_sea_surface(time, ts, cn, iceh)
type (time_type),                         intent(in)  :: time
real, dimension(:, :),                    intent(out) :: ts
real, dimension(size(ts,1),size(ts,2),2), intent(out) :: cn
real, dimension(size(ts,1),size(ts,2)),   intent(out) :: iceh

real, dimension(size(ts,1),size(ts,2))                :: sst, icec

real ::  t_sw_freeze = -1.8

!! TK Comment.  Deleted code associated with mcm_ice namelist 
!!              parameter that was in ice_spec version.
!!              It will not be available here...
!

  icec = 0.0; iceh = 0.0; sst = t_sw_freeze;
!  TK Mod: change these to ocean grid:
  call data_override('OCN', 'sic_obs', icec, time)
  call data_override('OCN', 'sit_obs', iceh, time)
!  LAT Mod: 3/6/09 - added override for restoring historical SST case
  if (do_restore_histsst) then
    call data_override('OCN', 'sst_histobs', sst, time)
  else
    call data_override('OCN', 'sst_obs', sst, time)
  endif
!  call data_override('ICE', 'sic_obs', icec, time)
!  call data_override('ICE', 'sit_obs', iceh, time)
!  call data_override('ICE', 'sst_obs', sst, time)

    where (icec >= 0.2)
      iceh = max(iceh, 1.0)
! TK Removed following statement.  It was not in earlier 
!       (pre-Inchon) code... 4/3/03
!!       icec = 1.0
      sst = t_sw_freeze
    else where
      icec = 0.0
      iceh = 0.0
    end where
  
! TK Mod: included following statement here so it is done
!   even for non mcm_ice...
  
  where (icec==0.0 .and. sst<=t_sw_freeze) sst = t_sw_freeze+1e-10

  cn(:,:,2) = icec
  cn(:,:,1) = 1-cn(:,:,2)
  ts = sst+Tfreeze

  return
end subroutine get_sea_surface

!#######################################################################


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! get_heat_flux_adj - get climatological heat flux adjustment from restoring  run           !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine get_heat_flux_adj(time, qflux_adj)

! TK NOTE: This routine outputs qflux_adj.  This is the climatological value of the heat
!          flux restoring for the given time in the seasonal cycle.  It is also interpolated
!          to the appropriate model grid.  The output field is qflux_adj
!          to avoid confusing it with hflx in the main ocean mixed layer program.
!          8/22/01.

! TK NOTE:  Modified to work with data override routines.
!          4/14/03

type (time_type),                         intent(in)  :: time
real, dimension(:, :),                    intent(out) :: qflux_adj

!  TK Mod: input data on ocean grid:
  call data_override('OCN', 'qflux_adj', qflux_adj, time)
  
  return
end subroutine get_heat_flux_adj

!#######################################################################
! dummy routine

subroutine ocean_stock_pe(Ocean_state, index, value, time_index)
! Right now this is a dummy subroutine.
use stock_constants_mod, only : ISTOCK_WATER, ISTOCK_HEAT
type(ocean_state_type), pointer     :: Ocean_state    ! internal ocean state
integer,               intent(in)   :: index
real,                  intent(out)  :: value
integer, optional,     intent(in)   :: time_index

  integer :: i, j, k

  value = 0.0

  if (.not.associated(Ocean_state)) return
  if (.not.Ocean_state%is_ocean_pe) return

  if(.not.stock_warning_issued) then
     call error_mesg('ocean_stock_pe','Stocks not yet implemented. Returning zero.',WARNING)
     stock_warning_issued = .true.
  endif

  select case (index)
    case (ISTOCK_WATER) ; value = 0.0
    case (ISTOCK_HEAT) ; value = 0.0
    case default ; value = 0.0
  end select

end subroutine ocean_stock_pe
!#######################################################################
subroutine ocean_model_data2D_get(OS,Ocean, name, array2D,isc,jsc)
  type(ocean_state_type),     pointer    :: OS
  type(ocean_public_type),    intent(in) :: Ocean
  character(len=*)          , intent(in) :: name
  real, dimension(isc:,jsc:), intent(out):: array2D
  integer                   , intent(in) :: isc,jsc
  
  array2D(isc:,jsc:) = 0.0
  
end subroutine ocean_model_data2D_get

subroutine ocean_model_data1D_get(OS,Ocean, name, value)
  type(ocean_state_type),     pointer    :: OS
  type(ocean_public_type),    intent(in) :: Ocean
  character(len=*)          , intent(in) :: name
  real                      , intent(out):: value

  value = 0.0

end subroutine ocean_model_data1D_get

! Subroutines for calculating checksums
SUBROUTINE ocean_public_type_chksum(id, timestep, bnd_type)
  USE fms_mod, ONLY: stdout
  USE mpp_mod, ONLY: mpp_chksum

  CHARACTER(len=*), INTENT(in) :: id
  INTEGER, INTENT(in) :: timestep
  TYPE(ocean_public_type), INTENT(in) :: bnd_type

  INTEGER :: n, m, outunit

  outunit = stdout()
100 FORMAT("CHECKSUM::",A32," = ",Z20)
101 FORMAT("CHECKSUM::",A16,A,"%",A," = ",Z20)

  WRITE (outunit,*) 'BEGIN CHECKSUM(ocean_public_type):: ', id, timestep
  WRITE (outunit,100) 'ocean_public_type%t_surf ', mpp_chksum(bnd_type%t_surf)
  WRITE (outunit,100) 'ocean_public_type%s_surf ', mpp_chksum(bnd_type%s_surf)
  WRITE (outunit,100) 'ocean_public_type%sea_lev', mpp_chksum(bnd_type%sea_lev)
  WRITE (outunit,100) 'ocean_public_type%frazil ', mpp_chksum(bnd_type%frazil)
  WRITE (outunit,100) 'ocean_public_type%u_surf ', mpp_chksum(bnd_type%u_surf)
  WRITE (outunit,100) 'ocean_public_type%area   ', mpp_chksum(bnd_type%area)
  WRITE (outunit,100) 'ocean_public_type%v_surf ', mpp_chksum(bnd_type%v_surf)

  DO n = 1, bnd_type%fields%num_bcs
     DO m = 1, bnd_type%fields%bc(n)%num_fields
        WRITE (outunit,101) 'ocean_public_type%', TRIM(bnd_type%fields%bc(n)%name),&
             & TRIM(bnd_type%fields%bc(n)%field(m)%name),&
             & mpp_chksum(bnd_type%fields%bc(n)%field(m)%values)
     END DO
  END DO
END SUBROUTINE ocean_public_type_chksum

SUBROUTINE ice_ocn_bnd_type_chksum(id, timestep, bnd_type)
  USE fms_mod, ONLY: stdout
  USE mpp_mod, ONLY: mpp_chksum

  CHARACTER(len=*), INTENT(in) :: id
  INTEGER, INTENT(in) :: timestep
  TYPE(ice_ocean_boundary_type), INTENT(in) :: bnd_type

  INTEGER :: n, m, outunit

  outunit = stdout()
100 FORMAT("CHECKSUM::",A32," = ",Z20)
101 FORMAT("CHECKSUM::",A16,A,"%",A," = ",Z20)

  WRITE (outunit,*) 'BEGIN CHECKSUM(ice_ocean_boundary_type):: ', id, timestep
  WRITE (outunit,100) 'ice_ocean_boundary_type%u_flux         ', mpp_chksum(bnd_type%u_flux)
  WRITE (outunit,100) 'ice_ocean_boundary_type%v_flux         ', mpp_chksum(bnd_type%v_flux)
  WRITE (outunit,100) 'ice_ocean_boundary_type%t_flux         ', mpp_chksum(bnd_type%t_flux)
  WRITE (outunit,100) 'ice_ocean_boundary_type%q_flux         ', mpp_chksum(bnd_type%q_flux)
  WRITE (outunit,100) 'ice_ocean_boundary_type%salt_flux      ', mpp_chksum(bnd_type%salt_flux)
  WRITE (outunit,100) 'ice_ocean_boundary_type%lw_flux        ', mpp_chksum(bnd_type%lw_flux)
  WRITE (outunit,100) 'ice_ocean_boundary_type%sw_flux_nir_dir', mpp_chksum(bnd_type%sw_flux_nir_dir)
  WRITE (outunit,100) 'ice_ocean_boundary_type%sw_flux_nir_dif', mpp_chksum(bnd_type%sw_flux_nir_dif)
  WRITE (outunit,100) 'ice_ocean_boundary_type%sw_flux_vis_dir', mpp_chksum(bnd_type%sw_flux_vis_dir)
  WRITE (outunit,100) 'ice_ocean_boundary_type%sw_flux_vis_dif', mpp_chksum(bnd_type%sw_flux_vis_dif)
  WRITE (outunit,100) 'ice_ocean_boundary_type%lprec          ', mpp_chksum(bnd_type%lprec)
  WRITE (outunit,100) 'ice_ocean_boundary_type%fprec          ', mpp_chksum(bnd_type%fprec)
  WRITE (outunit,100) 'ice_ocean_boundary_type%runoff         ', mpp_chksum(bnd_type%runoff)
  WRITE (outunit,100) 'ice_ocean_boundary_type%calving        ', mpp_chksum(bnd_type%calving)
  WRITE (outunit,100) 'ice_ocean_boundary_type%p              ', mpp_chksum(bnd_type%p)
  WRITE (outunit,100) 'ice_ocean_boundary_type%runoff_hflx    ', mpp_chksum(bnd_type%runoff_hflx)
  WRITE (outunit,100) 'ice_ocean_boundary_type%calving_hflx   ', mpp_chksum(bnd_type%runoff_hflx)

  DO n = 1, bnd_type%fluxes%num_bcs
     DO m = 1, bnd_type%fluxes%bc(n)%num_fields
        WRITE (outunit,101) 'ice_ocean_boundary_type%', TRIM(bnd_type%fluxes%bc(n)%name),&
             & TRIM(bnd_type%fluxes%bc(n)%field(m)%name),&
             & mpp_chksum(bnd_type%fluxes%bc(n)%field(m)%values)
     END DO
  END DO
END SUBROUTINE ice_ocn_bnd_type_chksum

end module ocean_model_mod
