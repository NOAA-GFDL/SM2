
      module m_tracname_mod
!-----------------------------------------------------------
!       ... List of advected and non-advected trace species, and
!           surface fluxes for the advected species.
!-----------------------------------------------------------

      use mo_grid_mod,   only : pcnst
      use chem_mods_mod, only : grpcnt

      implicit none

character(len=128), parameter :: version     = '$Id$'
character(len=128), parameter :: tagname     = '$Name$'
logical                       :: module_is_initialized = .false.

      save

      character(len=8) :: tracnam(pcnst)          ! species names
      character(len=8) :: natsnam(max(1,grpcnt))  ! names of non-advected trace species

      end module m_tracname_mod
