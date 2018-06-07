                    module sea_esf_rad_mod

! <CONTACT EMAIL="Fei.Liu@noaa.gov">
!   fil
! </CONTACT>
! <REVIEWER EMAIL="Stuart.Freidenreich@noaa.gov">
!   smf
! </REVIEWER>
! <REVIEWER EMAIL="">
!   
! </REVIEWER>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!   Code to initialize, commpute, and clean up radiation calculation. 
! </OVERVIEW>
! <DESCRIPTION>
!   The radiation component that initializes, deployes, and ends longwave,
!   shortwave, and diagnostics calculation in the FMS model.
! </DESCRIPTION>

!  shared modules:

use mpp_mod,              only: input_nml_file
use fms_mod,              only: open_namelist_file, fms_init, &
                                mpp_pe, mpp_root_pe, stdlog, &
                                file_exist, write_version_number, &
                                check_nml_error, error_mesg, &
                                FATAL, close_file, &
                                mpp_clock_id, mpp_clock_begin, &
                                mpp_clock_end, CLOCK_ROUTINE, &
                                CLOCK_MODULE
use time_manager_mod,     only: time_manager_init, time_type

!  shared radiation package modules:

use rad_utilities_mod,    only: radiation_control_type, &
                                radiative_gases_type, & 
                                astronomy_type

use aerosolrad_types_mod, only: aerosolrad_control_type

!   radiation package modules:

use longwave_driver_mod,  only: longwave_driver_init,   &
                                longwave_driver_time_vary, &
                                longwave_driver, &
                                longwave_driver_endts, &
                                longwave_driver_end, &
                                lw_output_type, &
                                lw_diagnostics_type, &
                                lw_table_type, &
                                longwave_number_of_bands, &
                                longwave_get_tables, &
                                longwave_output_alloc, &
                                longwave_diag_alloc, &
                                longwave_dealloc, &
                                assignment(=)

use shortwave_driver_mod, only: shortwave_driver_init,  &
                                shortwave_driver,  &
                                shortwave_driver_end, &
                                shortwave_driver_time_vary, &
                                shortwave_number_of_bands, &
                                get_solar_constant, &
                                sw_output_type, &
                                shortwave_output_alloc, &
                                shortwave_output_dealloc, &
                                assignment(=)
                                

!----------------------------------------------------------------------

implicit none 
private 

!-----------------------------------------------------------------------
!    sea_esf_rad_mod is the driver for the sea_esf_rad radiation 
!    package.
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!------------ version number for this module ---------------------------

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'


!--------------------------------------------------------------------
!-- interfaces -----

public       &
            sea_esf_rad_init, sea_esf_rad, sea_esf_rad_time_vary,  &
            sea_esf_rad_endts,  sea_esf_rad_end

! inherited from longwave & shortwave modules
public  longwave_number_of_bands, longwave_get_tables, &
        longwave_output_alloc, longwave_dealloc, &
        lw_diagnostics_type, longwave_diag_alloc, &
        shortwave_output_alloc, shortwave_output_dealloc, &
        shortwave_number_of_bands, get_solar_constant, &
        lw_table_type, &
        sw_output_type, lw_output_type, assignment(=)

!---------------------------------------------------------------------
!--- namelist ---

logical :: dummy


namelist /sea_esf_rad_nml/   &
                            dummy

!---------------------------------------------------------------------
!---- public data ----


!---------------------------------------------------------------------
!---- private data ----


logical :: module_is_initialized = .false.    ! module initialized ?
integer :: longwave_clock, shortwave_clock    ! timing clocks


!---------------------------------------------------------------------
!---------------------------------------------------------------------



                         contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!######################################################################
! <SUBROUTINE NAME="sea_esf_rad_init">
!   <OVERVIEW>
!     Routine to initialize the radiation calculation
!   </OVERVIEW>
!   <DESCRIPTION>
!     This routine initializes the utilities and radiation utilities
!     modules. Then it reads in the radiation namelist from the input
!     namelist file and log the namelist in an output log file.
!   </DESCRIPTION>
!   <TEMPLATE>
!     CALL sea_esf_rad_init (pref_r, Rad_control)
!   </TEMPLATE>
!
!   <IN NAME="pref_r" TYPE="real">
!     Array containing two reference pressure profiles 
!     on the radiation grid for use in defining 
!     transmission functions in [pascals]
!   </IN>
! </SUBROUTINE>
subroutine sea_esf_rad_init (pref_r, Rad_control)

!---------------------------------------------------------------------
!   sea_esf_rad_init is the constructor for sea_esf_rad_mod.
!---------------------------------------------------------------------
real,  dimension(:,:),        intent(in)    :: pref_r
type(radiation_control_type), intent(inout) :: Rad_control

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       pref_r    array containing two reference pressure profiles 
!                 on the radiation grid for use in defining 
!                 transmission functions 
!                 [pascals]
!
!----------------------------------------------------------------------

!-------------------------------------------------------------------
!  local variables

      integer                           :: unit, io, ierr, logunit

!---------------------------------------------------------------------
!  local variables:
!
!        unit            io unit number used for namelist file
!        ierr            error code
!        io              error status returned from io operation
!        end
!        Lw_tables
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call fms_init
      call time_manager_init

!-----------------------------------------------------------------------
!    read namelist.
!-----------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=sea_esf_rad_nml, iostat=io)
      ierr = check_nml_error(io,'sea_esf_rad_nml')
#else   
      if ( file_exist('input.nml')) then
        unit =  open_namelist_file ( )
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=sea_esf_rad_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'sea_esf_rad_nml')
        end do
10      call close_file (unit)
      endif
#endif

!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
                           write (logunit, nml=sea_esf_rad_nml)

!---------------------------------------------------------------------
!    initialize the modules called by this module.
!---------------------------------------------------------------------
      call longwave_driver_init  (pref_r, Rad_control)
      call shortwave_driver_init (Rad_control)

!---------------------------------------------------------------------
!    initialize clocks to time various modules called by this module.
!---------------------------------------------------------------------
      longwave_clock =      &
                  mpp_clock_id ('   Physics_down: Radiation: lw', &
                        grain=CLOCK_ROUTINE)
      shortwave_clock =     &
                  mpp_clock_id ('   Physics_down: Radiation: sw', &
                        grain=CLOCK_ROUTINE)

!-------------------------------------------------------------------
!    mark the module as initialized.
!-------------------------------------------------------------------
      module_is_initialized = .true.

!-------------------------------------------------------------------


end subroutine sea_esf_rad_init


!####################################################################
 
subroutine sea_esf_rad_time_vary (Time, Rad_gases_tv)


!----------------------------------------------------------------------
type(time_type), intent(in)  :: Time
type(radiative_gases_type), intent(inout) :: Rad_gases_tv

 
      call shortwave_driver_time_vary (Time)
      call longwave_driver_time_vary  (Rad_gases_tv)
 

end subroutine sea_esf_rad_time_vary
 

!#######################################################################        ######

subroutine sea_esf_rad_endts (Rad_gases_tv)
 
type(radiative_gases_type), intent(in) :: Rad_gases_tv
 
    call longwave_driver_endts (Rad_gases_tv)

end subroutine sea_esf_rad_endts 



!#####################################################################
! <SUBROUTINE NAME="sea_esf_rad">
!   
!   <OVERVIEW>
!     The radiation component interface of the climate model
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine calls longwave radiation computation subroutine, 
!     shortwave radiation computation subroutine, radiation diagnostics
!     computation routine, and finally it deallocates all previously
!     allocated memory spaces of temporary arrays.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call sea_esf_rad (Atmos_input, Surface, Astro, Rad_gases, &
!                       aerooptdep, aeroasymfac, aerosctopdep, aeroextopdep, &
!                       crndlw, cmxolw, emrndlw, emmxolw, camtsw, cldsct, cldext, cldasymm, &
!                       Lw_output, Sw_output)
!   </TEMPLATE>
!
!   <IN NAME="Atmos_input" TYPE="atmos_input_type">
!     Atmos_input_type variable containing the atmospheric
!     input fields on the radiation grid 
!   </IN>
!   <IN NAME="Astro" TYPE="astronomy_type">
!     Astronomy_type variable containing the astronomical
!     input fields on the radiation grid  
!   </IN>
!   <IN NAME="Rad_gases" TYPE="radiative_gases_type">
!     Radiative_gases_type variable containing the radiative 
!     gas input fields on the radiation grid 
!   </IN>
!   <IN NAME="aerooptdep" TYPE="real">
!    Longwave aerosol optical depth (mean in model layers)
!   </IN>
!   <IN NAME="crndlw" TYPE="real">
!    Longwave cloud amount for random overlap clouds
!   </IN>
!   <IN NAME="cmxolw" TYPE="real">
!    Longwave cloud amount for maximum overlap clouds
!   </IN>
!   <IN NAME="emrndlw" TYPE="real">
!    cloud emissivity for random overlap clouds
!    by longwave band and profile
!   </IN>
!   <IN NAME="emmxolw" TYPE="real">
!    cloud emissivity for maximum overlap clouds
!    by longwave band and profile
!   </IN>
!   <IN NAME="camtsw" TYPE="real">
!    Cloud amount for shortwave clouds. If stochastic clouds is implemented
!    then cloud amount by band.
!   </IN>
!   <IN NAME="cldext" TYPE="real">
!    Cloud extinction parameter
!   </IN>
!   <IN NAME="cldsct" TYPE="real">
!    Cloud single scattering albedo
!   </IN>
!   <IN NAME="cldasymm" TYPE="real">
!    Cloud asymmetric parameter
!   </IN>
!   <INOUT NAME="Lw_output" TYPE="lw_output_type">
!     The longwave radiation calculation result
!   </INOUT>
!   <INOUT NAME="Sw_output" TYPE="sw_output_type">
!     The shortwave radiation calculation result
!   </INOUT>
!   <IN NAME="Surface" TYPE="surface_type">
!    Surface data as boundary condition to radiation
!   </IN>
! </SUBROUTINE>

subroutine sea_esf_rad (press, pflux, temp, tflux, rh2o, deltaz, &
                        asfc_vis_dir, asfc_nir_dir, &
                        asfc_vis_dif, asfc_nir_dif, &
                        Astro, Rad_gases, aerooptdep, aerooptdep_volc, &
                        aeroasymfac, aerosctopdep, aeroextopdep, &
                        crndlw, cmxolw, emrndlw, emmxolw, &
                        camtsw, cldsct, cldext, cldasymm, &
                        flag_stoch, Rad_control, Aerosolrad_control, &
                        Lw_output, Sw_output, Lw_diagnostics)

!-----------------------------------------------------------------------
!     sea_esf_rad calls the modules which calculate the long- and short-
!     wave radiational heating terms and fluxes and the radiation diag-
!     nostics module which provides radiation package diagnostics.
!-----------------------------------------------------------------------

real, dimension(:,:,:),       intent(in)     :: press, pflux, temp, &
                                                tflux, rh2o, deltaz
real, dimension(:,:),         intent(in)     :: asfc_vis_dir, &
                                                asfc_nir_dir, &
                                                asfc_vis_dif, &
                                                asfc_nir_dif
type(astronomy_type),         intent(in)     :: Astro
type(radiative_gases_type),   intent(inout)  :: Rad_gases
real,dimension(:,:,:,:),      intent(in)     :: aerooptdep, &
                                                aerooptdep_volc
real,dimension(:,:,:,:),      intent(in)     :: aeroasymfac, &
                                                aerosctopdep, &
                                                aeroextopdep
real, dimension(:,:,:,:),     intent(in)     :: crndlw
real, dimension(:,:,:),       intent(in)     :: cmxolw
real, dimension(:,:,:,:,:),   intent(in)     :: emrndlw, emmxolw
real, dimension(:,:,:,:),     intent(in)     :: camtsw
real, dimension(:,:,:,:,:),   intent(in)     :: cldsct, cldext, cldasymm
integer,                      intent(in)     :: flag_stoch
type(radiation_control_type),  intent(in)    :: Rad_control
type(aerosolrad_control_type), intent(in)    :: Aerosolrad_control
type(lw_output_type), dimension(:), intent(inout)  :: Lw_output
type(sw_output_type), dimension(:), intent(inout)  :: Sw_output 
type(lw_diagnostics_type)         , intent(inout)  :: Lw_diagnostics
!---------------------------------------------------------------------
!  intent(in) variables:
!
!      Atmos_input   atmospheric input fields          
!                    [ atmos_input_type ]
!      Surface       surface variables 
!                    [ surface_type ]
!      Astro         astronomical input fields            
!                    [ astronomy_type ]
!      Rad_gases     radiative gas input fields   
!                    [ radiative_gases_type ]
!      aerooptdep    longwave aerosol optical depth in model layers
!      aeroasymfac, aerosctopdep, aeroextopdep
!                     shortwave aerosol radiative properties
!      crndlw       longwave cloud amount for random overlap clouds
!      cmxolw       longwave cloud amount for maximum overlap clouds
!      emrndlw      longwave cloud emissivity for random overlap clouds by band
!      emmxolw      longwave cloud emissivity for maximum overlap clouds by band
!      camtsw       shortwave cloud amount,
!                   if stochastic clouds is implemented then cloud amount by band
!      cldsct, cldext, cldasymm
!                   shortwave cloud radiative properties
!
!  intent(out) variables:
!
!      Lw_output     longwave radiation output data from the sea_esf_rad
!                    radiation package
!                    [ lw_output_type ]
!      Sw_output     shortwave radiation output data from the
!                    sea_esf_rad radiation package 
!                    [ sw_output_type ]
!      Lw_diagnostics    used to hold desired diagnostics from 
!                        longwave_driver_mod so they may be passed
!                        to radiation_diag_mod
!                        [ lw_diagnostics_type ]
!
!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('sea_esf_rad_mod',   &
              'module has not been initialized', FATAL )
      endif

!----------------------------------------------------------------------
!    compute longwave radiation.
!----------------------------------------------------------------------
    if (Rad_control%do_lw_rad) then
      call mpp_clock_begin (longwave_clock)
      call longwave_driver (press, pflux, temp, tflux, rh2o, deltaz,  &
                            Rad_gases, emrndlw, emmxolw, crndlw, cmxolw, &
                            aerooptdep, aerooptdep_volc, &
                            flag_stoch, Rad_control, Aerosolrad_control, &
                            Lw_output, Lw_diagnostics)
      call mpp_clock_end (longwave_clock)
    endif

!----------------------------------------------------------------------
!    compute shortwave radiation.
!----------------------------------------------------------------------
    if (Rad_control%do_sw_rad) then
      call mpp_clock_begin (shortwave_clock)
      call shortwave_driver (press, pflux, temp, rh2o, deltaz, &
                             asfc_vis_dir, asfc_nir_dir, &
                             asfc_vis_dif, asfc_nir_dif, Astro, &
                             aeroasymfac, aerosctopdep, aeroextopdep, &
                             Rad_gases, camtsw, cldsct, cldext, cldasymm, &
                             flag_stoch, Rad_control, Aerosolrad_control, Sw_output)
      call mpp_clock_end (shortwave_clock)
    endif

!--------------------------------------------------------------------

end subroutine sea_esf_rad


!###################################################################
! <SUBROUTINE NAME="sea_esf_rad_end">
! 
!   <OVERVIEW>
!     Ends radiation calculation.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine ends longwave, shortwave, and radiation
!     diagnostics calculation.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call sea_esf_rad_end
!   </TEMPLATE>
! </SUBROUTINE>

subroutine sea_esf_rad_end (Rad_control)
 
!-------------------------------------------------------------------
!    sea_esf_rad_end is the destructor for the sea_esf_rad module.
!-------------------------------------------------------------------
type(radiation_control_type), intent(in) :: Rad_control

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('sea_esf_rad_mod',   &
              'module has not been initialized', FATAL )
      endif

!--------------------------------------------------------------------
!    close out the modules initialized by this module.
!--------------------------------------------------------------------
      call longwave_driver_end
      call shortwave_driver_end (Rad_control)
 
!--------------------------------------------------------------------
!    mark the module as uninitialized.
!--------------------------------------------------------------------
      module_is_initialized = .false.

!--------------------------------------------------------------------


end subroutine sea_esf_rad_end


!####################################################################
      

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!####################################################################



                 end module sea_esf_rad_mod

