
module cloudrad_types_mod

!--------------------------------------------------------------------

implicit none 
private 

!--------------------------------------------------------------------
!---- public data structures ----

public cld_specification_type

type cld_specification_type
   real, dimension(:,:,:,:),  pointer :: tau=>NULL(),  &
                                         camtsw_band=>NULL(), &
                                         crndlw_band=>NULL(), &
                                         lwp_lw_band=>NULL(), &
                                         iwp_lw_band=>NULL(), &
                                         lwp_sw_band=>NULL(), &
                                         iwp_sw_band=>NULL(), &
                                         reff_liq_lw_band=>NULL(),   &
                                         reff_ice_lw_band=>NULL(), &
                                         reff_liq_sw_band=>NULL(),   &
                                         reff_ice_sw_band=>NULL()
   real, dimension(:,:,:),    pointer :: lwp=>NULL(),   &
                                         iwp=>NULL(),  &
                                         reff_liq=>NULL(),   &
                                         reff_ice=>NULL(), &
                                         reff_liq_lim=>NULL(),   &
                                         reff_ice_lim=>NULL(), &
                                         liq_frac=>NULL(), &
                                         cloud_water=>NULL(), &
                                         cloud_ice=>NULL(),  &
                                         cloud_area=>NULL(), &
                                         cloud_droplet=>NULL(), &
                                         cloud_ice_num=>NULL(), &
                                         rain =>NULL(), &
                                         snow =>NULL(), &
                                         rain_size =>NULL(), &
                                         snow_size =>NULL(), &
                                         reff_liq_micro=>NULL(),   &
                                         reff_ice_micro=>NULL(),&
                                         camtsw=>NULL(),   &
                                         cmxolw=>NULL(),  &
                                         crndlw=>NULL()
   integer, dimension(:,:,:), pointer :: cld_thickness=>NULL()
   integer, dimension(:,:,:,:), pointer :: stoch_cloud_type=>NULL()
   integer, dimension(:,:,:,:), pointer :: cld_thickness_lw_band=>NULL()
   integer, dimension(:,:,:,:), pointer :: cld_thickness_sw_band=>NULL()
   integer, dimension(:,:),   pointer :: ncldsw=>NULL(),   &
                                         nmxolw=>NULL(),&
                                         nrndlw=>NULL()
   integer, dimension(:,:,:), pointer :: ncldsw_band=>NULL(),   &
                                         nrndlw_band=>NULL()
   logical, dimension(:,:,:), pointer :: hi_cloud=>NULL(),   &
                                         mid_cloud=>NULL(),  &
                                         low_cloud=>NULL(),   &
                                         ice_cloud=>NULL()
end type cld_specification_type

!--------------------------------------------------------------------

public cldrad_properties_type

!    cldext     shortwave cloud extinction coefficient by band [km**(-1)]
!    cldasymm   shortwave cloud asymmetry factor by band [dimensionless]
!    cldsct     shortwave cloud scattering coefficient by band [km**(-1)]
!    emmxolw    longwave cloud emissivity for maximum overlap clouds by band [dimensionless]
!    emrndlw    longwave cloud emissivity for random overlap clouds by band [dimensionless]
!    abscoeff   longwave cloud absorption coefficient  by band [km**(-1)]
!  Following are no longer used:
!    cldemiss   
!    cirabsw
!    cirrfsw
!    cvisrfsw

type cldrad_properties_type
     real, dimension(:,:,:,:,:), pointer :: cldext=>NULL(),   &
                                            cldasymm=>NULL(), &
                                            cldsct=>NULL()
     real, dimension(:,:,:,:,:), pointer :: emmxolw=>NULL(),  &
                                            emrndlw=>NULL(),  &
                                            abscoeff=>NULL(), &
                                            cldemiss=>NULL()
     real, dimension(:,:,:),     pointer :: cirabsw=>NULL(), &
                                            cirrfsw=>NULL(), &
                                            cvisrfsw=>NULL()
end type cldrad_properties_type

!--------------------------------------------------------------------

public cloudrad_control_type

type cloudrad_control_type
    logical :: do_pred_cld_microphys
    logical :: do_presc_cld_microphys
    logical :: do_bulk_microphys
    logical :: do_sw_micro
    logical :: do_lw_micro
    logical :: do_strat_clouds
    logical :: do_no_clouds
    logical :: do_donner_deep_clouds
    logical :: do_uw_clouds
    logical :: do_random_overlap
    logical :: do_max_random_overlap
    logical :: do_stochastic_clouds
    logical :: use_temp_for_seed
    logical :: do_ica_calcs
    logical :: do_liq_num
    logical :: do_ice_num
    logical :: using_fu2007
    integer :: num_sw_cloud_bands
    integer :: num_lw_cloud_bands
end type cloudrad_control_type

!--------------------------------------------------------------------

public microphysics_type

type microphysics_type
character(len=64) :: scheme_name
real, dimension(:,:,:), pointer    :: conc_ice=>NULL(),   &
                                      conc_drop=>NULL(),      &
                                      size_ice=>NULL(),   &
                                      size_drop=>NULL(),     &
                                      size_snow=>NULL(),   &
                                      conc_snow=>NULL(),     &
                                      size_rain=>NULL(),     &
                                      conc_rain=>NULL(),   &
                                      cldamt=>NULL(),      &
                                      droplet_number=>NULL(), &
                                      ice_number=>NULL()
!  The following are activated for stochastic clouds
real, dimension(:,:,:,:), pointer :: stoch_conc_ice=>NULL(),   &
                                     stoch_conc_drop=>NULL(),  &
                                     stoch_size_ice=>NULL(),   &
                                     stoch_size_drop=>NULL(),  &
                                     stoch_cldamt=>NULL(),     &
                                     stoch_droplet_number=>NULL(), &
                                     stoch_ice_number=>NULL()
integer, dimension(:,:,:,:), pointer ::  stoch_cloud_type=>NULL()

!  In practice, we allocate a single set of columns for the
!  stochastic clouds, then point to sections of the larger array
!  with the lw_ and sw_ pointer arrays. 
!  i.e., lw_stoch_conc_ice => stoch_conc_ice(:, :, :, 1:numLwBands)

real, dimension(:,:,:,:), pointer :: lw_stoch_conc_ice=>NULL(),   &
                                     lw_stoch_conc_drop=>NULL(),  &
                                     lw_stoch_size_ice=>NULL(),   &
                                     lw_stoch_size_drop=>NULL(),  &
                                     lw_stoch_cldamt=>NULL(),     &
                                     lw_stoch_droplet_number=>NULL(), &
                                     lw_stoch_ice_number=>NULL(), &
                                     sw_stoch_conc_ice=>NULL(),   &
                                     sw_stoch_conc_drop=>NULL(),  &
                                     sw_stoch_size_ice=>NULL(),   &
                                     sw_stoch_size_drop=>NULL(),  &
                                     sw_stoch_cldamt=>NULL(),     &
                                     sw_stoch_droplet_number=>NULL(), &
                                     sw_stoch_ice_number=>NULL()

!  The following are only activated for large-scale cloud diagnostics

real, dimension(:,:,:), pointer    :: lsc_conc_drop=>NULL(),      &
                                      lsc_size_drop=>NULL(),     &
                                      lsc_cldamt=>NULL(),      &
                                      lsc_droplet_number=>NULL()
end type microphysics_type

!--------------------------------------------------------------------

public microrad_properties_type

!    cldext
!    cldsct
!    cldasymm
!    abscoeff

type microrad_properties_type
   character(len=64) :: scheme_name
   real, dimension(:,:,:,:), pointer :: cldext=>NULL(),  &
                                        cldsct=>NULL(), &
                                        cldasymm=>NULL(),    &
                                        abscoeff=>NULL()
end type microrad_properties_type

!####################################################################

end module cloudrad_types_mod

