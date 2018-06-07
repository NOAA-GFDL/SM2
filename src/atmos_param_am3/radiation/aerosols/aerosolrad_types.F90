 
module aerosolrad_types_mod

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!------ public data structures ------
!--------------------------------------------------------------------

public aerosolrad_control_type

type aerosolrad_control_type
    logical :: do_aerosol
    logical :: do_swaerosol
    logical :: do_lwaerosol
    logical :: volcanic_sw_aerosols
    logical :: volcanic_lw_aerosols
    logical :: do_swaerosol_forcing
    logical :: do_lwaerosol_forcing
    integer :: indx_swaf
    integer :: indx_lwaf
    integer :: num_sw_aerosol_bands
    integer :: num_lw_aerosol_bands
end type aerosolrad_control_type

!--------------------------------------------------------------------

public   aerosolrad_diag_type

type aerosolrad_diag_type
     real, dimension(:,:,:,:,:), pointer  :: extopdep=>NULL(), &
                                             absopdep=>NULL(), &
                                             asymdep=>NULL()

     real, dimension(:,:,:,:), pointer  :: extopdep_vlcno=>NULL(), &
                                           absopdep_vlcno=>NULL(), &
                                           sw_heating_vlcno=>NULL(), &
                                           lw_extopdep_vlcno=>NULL(), &
                                           lw_absopdep_vlcno=>NULL()

     real, dimension(:,:,:,:),   pointer  :: sw_ext=>NULL(), &
                                             sw_ssa=>NULL(), &
                                             sw_asy=>NULL(), &
                                             lw_ext=>NULL(), &
                                             lw_ssa=>NULL(), &
                                             lw_asy=>NULL()
end type aerosolrad_diag_type

!####################################################################

end module aerosolrad_types_mod

