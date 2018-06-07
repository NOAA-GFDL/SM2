 
module aerosol_types_mod

!--------------------------------------------------------------------

use time_manager_mod,  only: time_type
use interpolator_mod,  only: interpolate_type

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!------ public data structures ------
!--------------------------------------------------------------------

public   aerosol_type

type aerosol_type
     real,       dimension(:,:,:,:), pointer :: aerosol=>NULL()
     logical,    dimension(:,:),     pointer :: family_members=>NULL()
     character(len=64), dimension(:), pointer :: aerosol_names=>NULL()
end type aerosol_type

!--------------------------------------------------------------------

public aerosol_time_vary_type

!  Interp                interpolate_type variable containing the
!                        information about the aerosol species
!  Time                  time for which data is obtained from
!                        aerosol timeseries
!  being_overridden      is a given aerosol field to be overridden
!                        based on the model data_table?
!  output_override_info  should override info about each
!                        aerosol field be output (will be
!                        set to .false. after first time step)
!  override_counter      used to count calls to aerosol_endts
!                        so that output_override_info may
!                        set to .false. after physics_up
!  nfields               number of active aerosol species
!  nfamilies             number of active aerosol families

type aerosol_time_vary_type
     type(interpolate_type), dimension(:), pointer :: Interp=>NULL()
     type(time_type),        dimension(:), pointer :: Time=>NULL()
     logical, dimension(:), pointer :: being_overridden=>NULL()
     integer  :: nfields=0
     integer  :: nfamilies=0
     logical :: output_override_info = .true.
     integer :: override_counter = 0
     logical :: variable_is_initialized = .false.
end type aerosol_time_vary_type

!------------------------------------------------------------------

end module aerosol_types_mod

