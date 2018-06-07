 
               module rad_utilities_mod
!
! <CONTACT EMAIL="Fei.Liu@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="Stuart.Freidenreich@noaa.gov">
!  smf
! </REVIEWER>
! 
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!  Code to define the derived data types and provide table search routines.
! </OVERVIEW>
! <DESCRIPTION>
!  This code is used in the radiation code package as a helper module.
!  It defines many derived data types used in radiation calculation.
!  This code also provides table search routines and simple arithmatic
!  routines.
! </DESCRIPTION>
!
use fms_mod,            only : fms_init, &
                               mpp_pe, mpp_root_pe, stdlog, &
                               write_version_number, &
                               error_mesg, FATAL, close_file
use time_manager_mod,   only : time_type

!--------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!    rad_utilities_mod contains radiation table search routines,
!    some band averaging routines, and the derived-type variables 
!    used in the radiation package.
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

character(len=128)  :: version =  '$Id$'
character(len=128)  :: tagname =  '$Name$'

!---------------------------------------------------------------------
!-------  interfaces --------

public    rad_utilities_init, rad_utilities_end

!------------------------------------------------------------------

public astronomy_type
    
!    solar
!    cosz
!    fracday
!    rrsun

type astronomy_type
     real, dimension(:,:), pointer  :: solar=>NULL(),   &
                                       cosz=>NULL(),  &
                                       fracday=>NULL()
     real, dimension(:,:,:), pointer  :: solar_p=>NULL(),   &
                                       cosz_p=>NULL(),  &
                                       fracday_p=>NULL()
     real    :: rrsun
end type astronomy_type

!--------------------------------------------------------------------
 
public astronomy_inp_type
 
!    zenith_angle   specified zenith angles [ degrees ]
!    fracday        specified daylight fraction [ fraction ]
!    rrsun          specified earth-sun distance, normalized by mean
!                   distance 

type astronomy_inp_type
     real, dimension(:,:), pointer  :: zenith_angle=>NULL()
     real, dimension(:,:), pointer  :: fracday=>NULL()
     real                           :: rrsun
end type astronomy_inp_type

!--------------------------------------------------------------------

public atmos_input_type

!    press
!    temp
!    rh2o
!    zfull
!    pflux
!    tflux
!    deltaz
!    phalf
!    relhum
!    cloudtemp
!    clouddeltaz
!    cloudvapor
!    aerosoltemp
!    aerosolvapor
!    aerosolpress
!    aerosolrelhum
!    tracer_co2
!    g_rrvco2
!    tsfc
!    psfc

type atmos_input_type
     real, dimension(:,:,:), pointer :: press=>NULL(),   &
                                        temp=>NULL(), &
                                        rh2o=>NULL(),  &
                                        zfull=>NULL(),  &
                                        pflux=>NULL(), &
                                        tflux=>NULL(),  &
                                        deltaz=>NULL(),  &
                                        phalf=>NULL(),   &
                                        relhum=>NULL(), &
                                        cloudtemp=>NULL(),   &
                                        clouddeltaz=>NULL(), &
                                        cloudvapor=>NULL(), &
                                        aerosoltemp=>NULL(), &
                                        aerosolvapor=>NULL(), &
                                        aerosolpress=>NULL(), &
                                        aerosolrelhum=>NULL(), &
                                        tracer_co2 => NULL()
     real, dimension(:,:),   pointer :: tsfc=>NULL(),   &
                                        psfc=>NULL()              
     real                            :: g_rrvco2
end type atmos_input_type

!-------------------------------------------------------------------

!public cld_specification_type

!    tau
!    lwp
!    iwp
!    reff_liq
!    reff_ice
!    liq_frac
!    cloud_water
!    cloud_ice
!    cloud_area
!    cloud_droplet
!    reff_liq_micro
!    reff_ice_micro
!    camtsw
!    cmxolw
!    crndlw
!    cld_thickness
!    ncldsw
!    nmxolw
!    nrndlw
!    hi_cloud
!    mid_cloud
!    low_cloud
!    ice_cloud

!------------------------------------------------------------------

public radiation_control_type

type radiation_control_type
    logical  :: do_totcld_forcing
    logical  :: using_restart_file
    logical  :: do_sw_rad
    logical  :: do_lw_rad
    logical  :: renormalize_sw_fluxes
    logical  :: hires_coszen
    integer  :: nzens
  ! from longwave_control_type
    logical  :: do_h2o
    logical  :: do_o3 
    logical  :: do_ch4_lw
    logical  :: do_n2o_lw
    logical  :: do_cfc_lw
    logical  :: do_co2_lw
    logical  :: do_co2_10um
  ! from shortwave_control_type
    logical  :: do_esfsw
    logical  :: do_diurnal
    logical  :: do_annual
    logical  :: do_daily_mean
    logical  :: do_cmip_sw_diagnostics
end type radiation_control_type

!------------------------------------------------------------------

public   radiative_gases_type
 
!    qo3
!    rrvch4
!    rrvn2o
!    rrvco2
!    rrvf11
!    rrvf12
!    rrvf113
!    rrvf22
!    rf11air
!    rf12air
!    rf113air
!    rf22air
!    time_varying_co2
!    time_varying_f11
!    time_varying_f12
!    time_varying_f113
!    time_varying_f22
!    time_varying_ch4
!    time_varying_n2o

type radiative_gases_type
     real, dimension(:,:,:), pointer :: qo3=>NULL()
     real                            :: rrvch4, rrvn2o, rrvco2,    &
                                        rrvf11, rrvf12, rrvf113,  &
                                        rrvf22, rf11air, rf12air,  &
                                        rf113air, rf22air, &
                                        co2_for_last_tf_calc,  &
                                        co2_tf_offset, &
                                        co2_for_next_tf_calc, &
                                        ch4_for_last_tf_calc,  &
                                        ch4_tf_offset, &
                                        ch4_for_next_tf_calc, &
                                        n2o_for_last_tf_calc,  &
                                        n2o_tf_offset, &
                                        n2o_for_next_tf_calc, &
                                        ch4_for_tf_calc, &
                                        n2o_for_tf_calc, &
                                        co2_for_tf_calc
     logical                         :: time_varying_co2,  &
                                        time_varying_f11, &
                                        time_varying_f12,  &
                                        time_varying_f113, &
                                        time_varying_f22,  &
                                        time_varying_ch4, &
                                        time_varying_n2o, &
                                        use_ch4_for_tf_calc, &
                                        use_n2o_for_tf_calc, &
                                        use_co2_for_tf_calc, &
                                        use_model_supplied_co2
     type(time_type)                 :: Co2_time, Ch4_time, N2o_time
end type radiative_gases_type

!------------------------------------------------------------------

public rad_output_type

!    tdt_rad
!    tdt_rad_clr
!    tdtsw
!    tdtsw_clr
!    tdtlw
!    flux_sw_surf
!    flux_lw_surf
!    coszen_angle

type rad_output_type
     real, dimension(:,:,:,:), pointer :: tdt_rad=>NULL(),  &
                                        ufsw=>NULL(),  &
                                        dfsw=>NULL(),  &
                                        tdtsw=>NULL()  
     real, dimension(:,:,:,:), pointer :: tdt_rad_clr=>NULL(), &
                                        ufsw_clr=>NULL(),  &
                                        dfsw_clr=>NULL(),  &
                                        tdtsw_clr=>NULL()
                                        
     real, dimension(:,:,:), pointer :: tdtlw=>NULL()
     real, dimension(:,:,:), pointer :: flxnet=>NULL()
     real, dimension(:,:,:), pointer :: flxnetcf=>NULL()
     real, dimension(:,:,:), pointer :: tdtlw_clr=>NULL()
     real, dimension(:,:,:),   pointer :: flux_sw_surf=>NULL(), &
                                        flux_sw_surf_refl_dir=>NULL(), &
                                        flux_sw_surf_dir=>NULL(), &
                                        flux_sw_surf_dif=>NULL(), &
                                        flux_sw_down_vis_dir=>NULL(), &
                                        flux_sw_down_vis_dif=>NULL(), &
                                       flux_sw_down_total_dir=>NULL(), &
                                       flux_sw_down_total_dif=>NULL(), &
                                        flux_sw_vis=>NULL(), &
                                        flux_sw_vis_dir=>NULL(), &
                                        flux_sw_refl_vis_dir=>NULL(), &
                                        flux_sw_vis_dif=>NULL()
     real, dimension(:,:,:),   pointer :: flux_sw_down_vis_clr=>NULL(), &
                                  flux_sw_down_total_dir_clr=>NULL(), &
                                  flux_sw_down_total_dif_clr=>NULL()
     real, dimension(:,:),   pointer :: flux_lw_surf=>NULL(), &
                                        coszen_angle=>NULL()
     real, dimension(:,:,:), pointer :: extinction=>NULL()
end type rad_output_type

!-------------------------------------------------------------------

public surface_type

!    asfc
!    land

type surface_type
    real, dimension(:,:),   pointer ::  asfc=>NULL(),   &
                                        land=>NULL(),  &
                                        asfc_vis_dir=>NULL(), &
                                        asfc_nir_dir=>NULL(), &
                                        asfc_vis_dif=>NULL(), &
                                        asfc_nir_dif=>NULL()
end type surface_type
 
!-------------------------------------------------------------------
!--------------------------------------------------------------------
!------- private data ------


logical :: module_is_initialized=.false.   ! module is initialized ?


!---------------------------------------------------------------------
!---------------------------------------------------------------------


                           contains


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! <SUBROUTINE NAME="rad_utilities_init">
!  <OVERVIEW>
!   Subroutine to initialize radiation utility package.
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine initializes rad_utilities_nml.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call rad_utilities_init
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine rad_utilities_init

!---------------------------------------------------------------------
!    rad_utilities_init is the constructor for rad_uti;lities_mod.
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

!-----------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)

!-------------------------------------------------------------------
!    mark the module as initialized.
!-------------------------------------------------------------------
      module_is_initialized = .true.

!------------------------------------------------------------------

end subroutine rad_utilities_init


!#####################################################################

! <SUBROUTINE NAME="rad_utilities_end">
!  <OVERVIEW>
!   Subroutine to close out the radiation utility package.
!  </OVERVIEW>
!  <DESCRIPTION>
!   This subroutine is the destructor for rad_utilies_mod. it marks
!   the module as uninitialized.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call rad_utilities_end
!  </TEMPLATE>
! </SUBROUTINE>
!
subroutine rad_utilities_end

!--------------------------------------------------------------------
!    rad_utilites_end is the destructor for rad_utilities_mod.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
       module_is_initialized = .false.

!-------------------------------------------------------------------

end subroutine rad_utilities_end

!#####################################################################


                     end module rad_utilities_mod


