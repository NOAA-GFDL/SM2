module atmosphere_mod
#include <fms_platform.h>

!-----------------------------------------------------------------------
!
!         interface for fv dynamical core and physics
!
!-----------------------------------------------------------------------
!---------------- m o d u l e   i n f o r m a t i o n ------------------


use time_manager_mod,    only: time_type, get_time, set_time, operator(+)

use fms_mod,             only: error_mesg, FATAL, stdlog,      &
                               write_version_number,           &
                               mpp_pe, mpp_root_pe,            &
                               set_domain,                     &
                               mpp_clock_id, mpp_clock_begin,  &
                               mpp_clock_end, CLOCK_SUBCOMPONENT, &
                               clock_flag_default


use    mpp_domains_mod, only: domain2d
use  field_manager_mod, only: MODEL_ATMOS
use      constants_mod, only: omega, cp_air, rdgas, grav, rvgas, kappa, radius, pstd_mks

use            fv_pack, only: nlon, mlat, nlev, beglat, endlat, beglon, &
                              endlon, rlonb, rlatb,  cold_start, ncnst, &
                              pnats, consv_te, ptop, fv_init, fv_domain, &
                              fv_end, change_time, restart_format, area, &
                              ak, bk, rlon, rlat, nt_prog, get_eta_level, &
                              esm2_bugs
use     fv_diagnostics, only: fv_diag_init, fv_diag, fv_time
use       timingModule, only: timing_on, timing_off
use  atmos_nudge_mod,   only: atmos_nudge_init, atmos_nudge_end
use fv_restart_mod,     only: fv_restart, write_fv_rst
use fv_dynamics_mod,    only: fv_dynamics
use fv_arrays_mod,      only: fv_print_chksums
use mpp_mod,            only: mpp_error, FATAL, NOTE
use tracer_manager_mod, only: get_number_tracers, get_tracer_index, NO_TRACER
use xgrid_mod,          only: grid_box_type

use update_fv_phys_mod, only: update_fv_phys
use block_control_mod,  only: block_control_type
use physics_driver_mod, only: surf_diff_type
use physics_types_mod,  only: physics_type, &
                              physics_tendency_type
use radiation_types_mod,only: radiation_type, compute_g_avg


!-----------------------------------------------------------------------

implicit none
private

!--- driver routines
public :: atmosphere_init, atmosphere_end, atmosphere_restart, &
          atmosphere_dynamics, atmosphere_state_update

!--- utility routines
public :: atmosphere_resolution, atmosphere_boundary, &
          atmosphere_grid_center, atmosphere_domain, &
          atmosphere_cell_area, atmosphere_control_data, &
          atmosphere_pref, &
          get_atmosphere_axes, get_bottom_mass, &
          get_bottom_wind, get_stock_pe, &
          set_atmosphere_pelist, reset_atmos_tracers

!--- physics/radiation data exchange routines
public :: atmos_radiation_driver_inputs, atmos_physics_driver_inputs

public  surf_diff_type

integer sec
integer seconds, days
integer liq_wat, ice_wat
!-----------------------------------------------------------------------

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'

!-----------------------------------------------------------------------
!---- namelist (saved in file input.nml) ----
!
! physics_window  The number of "i" by "j" rows processed each time
!                 the modular physics is called. To process the entire
!                 domain use physics_window = (/0,0/).
!                   [integer, default: physics_window = 0,0]

   integer, dimension(2) :: physics_window = (/0,0/)

! Note: the default size of window(1) is chosen to work for N30, N45, ...

!-----------------------------------------------------------------------
!---- private data ----

  type    (time_type) :: Time_step_atmos
  real                :: dt_atmos
  integer, dimension(4)              :: atmos_axes
  integer :: id_dynam, id_phys_down, id_phys_up, id_fv_diag

!-----------------------------------------------------------------------
  real, allocatable :: pref(:,:), dum1d(:)
  logical :: do_atmos_nudge


contains

 subroutine atmosphere_init (Time_init, Time, Time_step, Surf_diff, Grid_box)

 type (time_type),     intent(in)    :: Time_init, Time, Time_step
 type(surf_diff_type), intent(inout) :: Surf_diff
 type(grid_box_type),  intent(inout) :: Grid_box

  integer :: ss, ds

   call write_version_number ( version, tagname )

!---- compute physics/atmos time step in seconds ----

   Time_step_atmos = Time_step
   call get_time (Time_step_atmos, sec)
   dt_atmos = real(sec)

!----- initialize FV dynamical core -----

   call fv_init( sec )
   call fv_restart( days, seconds )

    if ( cold_start .or. change_time ) then
        fv_time = time
    else
        fv_time = set_time (seconds, days)
        call get_time (Time, ss,  ds)

        if( seconds /= ss .or. days /= ds )   call  error_mesg         &
            ('FV_init:','Time inconsistent between fv_rst and INPUT/atmos_model.res', FATAL)
    endif

!----- initialize atmos_axes and fv_dynamics diagnostics

    call fv_diag_init( atmos_axes, Time )
!    call fv_print_chksums( 'after fv_diag_init' )
   

!----- initialize physics interface -----
!----- initialize domains for reading global physics data -----

    call set_domain ( fv_domain )

!  --- initialize clocks for dynamics, physics_down and physics_up

    id_dynam     = mpp_clock_id ('FV dynamical core',   &
                       flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
    id_fv_diag   = mpp_clock_id ('FV Diag',   &
                       flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )

!--- allocate pref
    allocate(pref(nlev+1,2), dum1d(nlev+1))
!---------- reference profile -----------
    pref(nlev+1,1) = 101325.
    pref(nlev+1,2) = 81060.

    call get_eta_level ( nlev, pref(nlev+1,1), pref(1,1), dum1d )
    call get_eta_level ( nlev, pref(nlev+1,2), pref(1,2), dum1d )

!--- setup Grid_box area for physics
    allocate(Grid_box%area  (beglon:endlon, beglat:endlat))
    Grid_box%area  (beglon:endlon, beglat:endlat) = area (beglon:endlon, beglat:endlat)

!--- initialize nudging module ---
    call atmos_nudge_init ( Time, atmos_axes(1:3), flag=do_atmos_nudge )

    call fv_print_chksums( 'Exiting  atmosphere_init' )

    liq_wat = get_tracer_index(MODEL_ATMOS, 'liq_wat')
    ice_wat = get_tracer_index(MODEL_ATMOS, 'ice_wat')
 end subroutine atmosphere_init

 subroutine atmosphere_dynamics (Time,surf_diff)
#include "fv_arrays.h"
!
!        Time = time at the current time level
!        surf_diff = not used in subroutine; only included for compatability with AM3/AM4 switch
!
   type(time_type),intent(in)    :: Time
   type(surf_diff_type),intent(in), optional :: surf_diff !! This is included for capatability 
                                                          !!  and is not used in this subroutine.

   type(time_type) :: Time_prev, Time_next
   real zvir
#include "fv_point.inc"


   zvir = rvgas/rdgas - 1.
   call fv_print_chksums( 'Entering  dynamics' )

   Time_prev = Time                       ! two time-level scheme
   Time_next = Time + Time_step_atmos

!---- dynamics -----

   call timing_on('fv_dynamics')
   call mpp_clock_begin (id_dynam)
   call fv_dynamics (nlon,    mlat,   nlev,    beglat,   endlat,    &
        ncnst,   pnats,  .false., consv_te,            &
        u,       v,      delp,    pt,       q,         &
        ps,      pe,     pk,      pkz,      phis,      &
        omga,    peln,   ptop,    omega,    sec,       &
        zvir,    cp_air, rdgas,   kappa,  radius, ua, va, Time_next )
   call mpp_clock_end (id_dynam)
   call timing_off('fv_dynamics')

   call fv_print_chksums( 'Exiting  atmosphere_dynamics' )
 end subroutine atmosphere_dynamics


 subroutine atmosphere_end (Time, Grid_box)

 type (time_type), intent(in) :: Time
 type(grid_box_type),  intent(inout) :: Grid_box

!----- initialize domains for writing global physics data -----

    call set_domain ( fv_domain )
    call get_time (Time, seconds,  days)
    call write_fv_rst( 'RESTART/fv_rst.res', days, seconds, grav, &
         restart_format )

    call fv_end(days, seconds)      

    call atmos_nudge_end

    deallocate (pref, dum1d)

 end subroutine atmosphere_end

  !#######################################################################
  ! <SUBROUTINE NAME="atmosphere_restart">
  ! <DESCRIPTION>
  !  dummy routine.
  ! </DESCRIPTION>
  subroutine atmosphere_restart(timestamp)
    character(len=*),  intent(in) :: timestamp

    call error_mesg ('atmosphere_restart in atmosphere_mod', &
                     'intermediate restart capability is not implemented for this model', FATAL)

  end subroutine atmosphere_restart
  ! </SUBROUTINE>

!---------------------------------------------------------------
!    returns the number of longitude and latitude grid points
!    for either the local PEs grid (default) or the global grid
!---------------------------------------------------------------

 subroutine atmosphere_resolution (i_size, j_size, global)

   integer, intent(out)          :: i_size, j_size
   logical, intent(in), optional :: global
   logical :: local

   local = .TRUE.
   if( PRESENT(global) )local = .NOT.global

   if( local )then
       i_size = endlon - beglon + 1
       j_size = endlat - beglat + 1
   else
       i_size = nlon
       j_size = mlat
   end if
 end subroutine atmosphere_resolution


 subroutine atmosphere_pref (p_ref)
   real, dimension(:,:), intent(inout) :: p_ref

   p_ref = pref

 end subroutine atmosphere_pref

 subroutine atmosphere_control_data (i1, i2, j1, j2, kt, p_hydro, hydro, do_uni_zfull)
   integer, intent(out)           :: i1, i2, j1, j2, kt
   logical, intent(out), optional :: p_hydro, hydro, do_uni_zfull
   i1 = beglon
   i2 = endlon
   j1 = beglat
   j2 = endlat
   kt = nlev

!---non-hydrostatic is not supported by the fv-latlon core
   if (present(p_hydro)) p_hydro = .true.
   if (present(  hydro))   hydro = .true.

 end subroutine atmosphere_control_data

 subroutine atmosphere_cell_area  (area_out)
    real, dimension(:,:),  intent(out) :: area_out

    area_out(1:size(area_out,1), 1:size(area_out,2)) =  &
                                   area (beglon:endlon, beglat:endlat)

 end subroutine atmosphere_cell_area


!---this is a dummy routine needed for compatibility with the 
!---decoupling of physics and radiation from the dynamic cores
 subroutine atmosphere_grid_center (lon, lat)
!---------------------------------------------------------------
!    returns the longitude and latitude cell centers
!---------------------------------------------------------------
    real,    intent(out) :: lon(:,:), lat(:,:)   ! Unit: radian
! Local data:
    integer i,j

    do j = beglat,endlat
      do i = beglon,endlon
        lon(i-beglon+1,j-beglat+1) = rlon(i,j)
        lat(i-beglon+1,j-beglat+1) = rlat(i,j)
      end do
    end do

 end subroutine atmosphere_grid_center


!---------------------------------------------------------------
!    returns the longitude and latitude grid box edges
!    for either the local PEs grid (default) or the global grid
!---------------------------------------------------------------

 subroutine atmosphere_boundary (blon, blat, global)

    real,    intent(out)          :: blon(:,:), blat(:,:)   ! radian
    logical, intent(in), optional :: global
! Local:
    integer i,j
    logical :: local

    local = .TRUE.
    if( PRESENT(global) )local = .NOT.global

    if( local )then
        do i = beglon,endlon+1
           blon(i-beglon+1,:) = rlonb(i)
        end do
        do j = beglat,endlat+1
           blat(:,j-beglat+1) = rlatb(j)
        end do
    else
        do i=1,nlon+1
           blon(i,:) = rlonb(i)
        end do
        do j=1,mlat+1
           blat(:,j) = rlatb(j)
        end do
    end if

 end subroutine atmosphere_boundary


 subroutine set_atmosphere_pelist ()
!--- no-op for fv_latlon dy-core
 end subroutine set_atmosphere_pelist


 subroutine atmosphere_domain (Domain)
 type(domain2d), intent(inout) :: Domain

!  returns the domain2d variable associated with the coupling grid
!  note: coupling is done using the mass/temperature grid with no halos

   Domain = fv_domain

 end subroutine atmosphere_domain


 subroutine get_atmosphere_axes ( axes )
! returns the axis indices associated with the coupling grid

   integer, intent(out) :: axes (:)

!----- returns the axis indices for the atmospheric (mass) grid -----

     if ( size(axes(:)) < 0 .or. size(axes(:)) > 4 ) call error_mesg (    &
                           'get_atmosphere_axes in atmosphere_mod', &
                           'size of argument is incorrect', FATAL   )

     axes (1:size(axes(:))) = atmos_axes (1:size(axes(:)))

 
 end subroutine get_atmosphere_axes



 subroutine get_bottom_mass (t_bot, tr_bot, p_bot, z_bot, p_surf, slp)

#include "fv_arrays.h"
! returns temp, sphum, pres, height at the lowest model level
!         and surface pressure and sea level pressure

   real, intent(out), dimension(beglon:endlon,beglat:endlat)  &
        :: t_bot, p_bot, z_bot, p_surf, slp
   real, intent(out), dimension(beglon:endlon,beglat:endlat,ncnst-pnats):: tr_bot
   integer :: i, j, k, kr
   real zvir, rrg, sigtop, sigbot
   real, dimension(beglon:endlon,beglat:endlat) :: tref
   real, parameter :: tlaps = 6.5e-3
#include "fv_point.inc"

   rrg  = rdgas / grav
   zvir = rvgas/rdgas - 1.


   ! determine 0.8 sigma reference level
   sigtop = ak(1)/pstd_mks+bk(1)
   do k = 1, nlev 
      sigbot = ak(k+1)/pstd_mks+bk(k+1)
      if (sigbot+sigtop > 1.6) then
         kr = k
         exit
      endif   
      sigtop = sigbot
   enddo

!$omp parallel do default(shared)    &
!$omp private (i, j)
     do j = beglat, endlat
        do i = beglon, endlon
           p_surf(i,j) =  ps(i,j)
           t_bot(i,j) =  pt(i,j,nlev)

           p_bot(i,j) = delp(i,j,nlev)/(peln(i,nlev+1,j)-peln(i,nlev,j))
           z_bot(i,j) = rrg*t_bot(i,j)*(1.+zvir*q(i,j,nlev,1))*  &
                  (1. - pe(i,nlev,j)/p_bot(i,j))
           ! sea level pressure
           tref(i,j) = pt(i,j,kr)*(delp(i,j,kr)/((peln(i,kr+1,j)-peln(i,kr,j))*ps(i,j)))**(-rrg*tlaps)
           slp(i,j) = ps(i,j)*(1.+tlaps*phis(i,j)/(tref(i,j)*grav))**(1./(rrg*tlaps))
        enddo
     enddo
! Copy tracers
!$omp parallel do default(shared)    &
!$omp private (i, j, k)
     do k = 1,ncnst-pnats
        do j = beglat,endlat
           do i = beglon,endlon
              tr_bot(i,j,k) = q(i,j,nlev,k)
           enddo
        enddo
     enddo

 end subroutine get_bottom_mass



 subroutine get_bottom_wind (u_bot, v_bot)
#include "fv_arrays.h"
!-----------------------------------------------------------
! returns u and v on the mass grid at the lowest model level
!-----------------------------------------------------------

   real, intent(out), dimension(beglon:,beglat:) :: u_bot, v_bot

   integer i, j
#include "fv_point.inc"
!Balaji: this cannot work unless beglon=1: corrected declaration Lbounds
   do j=beglat,endlat
      do i=beglon,endlon
         u_bot(i,j) = u_srf(i,j)
         v_bot(i,j) = v_srf(i,j)
      enddo
   enddo

 end subroutine get_bottom_wind

 subroutine get_stock_pe(index, value)
#include "fv_arrays.h"
    integer, intent(in) :: index
    real, intent(out)   :: value
#ifdef USE_STOCK
    include 'stock.inc' 
#endif
    real wm(beglon:endlon, beglat:endlat)
    integer i,j,k
#include "fv_point.inc"
   
    select case (index)
#ifdef USE_STOCK
    case (ISTOCK_WATER)
#else
    case (1)
#endif
     
!----------------------
! Perform vertical sum:
!----------------------
     wm = 0.
     do j = beglat, endlat
        do k=1,nlev
           do i = beglon, endlon
! Note: There is a check in fv_pack.F90 that ensures that tracer number one
!       is sphum and that the cloud water and ice tracers exist.
              wm(i,j) = wm(i,j) + delp(i,j,k)*(q(i,j,k,1)+q(i,j,k,liq_wat)+q(i,j,k,ice_wat))
           enddo
        enddo
     enddo

!----------------------
! Horizontal sum:
!----------------------
     value = 0.
     do j = beglat, endlat
        do i = beglon, endlon
           value = value + wm(i,j)*area(i,j)
        enddo
     enddo
     value = value/grav

    case default
     value = 0.0
    end select

 end subroutine get_stock_pe 

 subroutine atmosphere_state_update (Time, Physics_tendency, Physics, Atm_block)
#include "fv_arrays.h"
   type(time_type),              intent(in) :: Time
   type (physics_tendency_type), intent(in) :: Physics_tendency
   type (physics_type),          intent(in) :: Physics
   type (block_control_type),    intent(in) :: Atm_block
   type(time_type) :: Time_next
!--- local variables ---
   integer:: ibs, ibe, jbs, jbe, nb
   real:: zvir
#include "fv_point.inc"

   zvir = rvgas/rdgas - 1.

!--- put u/v tendencies into haloed arrays u_dt and v_dt
!$OMP parallel do default(shared) private(nb, ibs, ibe, jbs, jbe)
   do nb = 1,Atm_block%nblks
     ibs = Atm_block%ibs(nb)
     ibe = Atm_block%ibe(nb)
     jbs = Atm_block%jbs(nb)
     jbe = Atm_block%jbe(nb)

     u_dt(ibs:ibe,jbs:jbe,:)   = Physics_tendency%block(nb)%u_dt
     v_dt(ibs:ibe,jbs:jbe,:)   = Physics_tendency%block(nb)%v_dt
     t_dt(ibs:ibe,jbs:jbe,:)   = Physics_tendency%block(nb)%t_dt
     q_dt(ibs:ibe,jbs:jbe,:,:) = Physics_tendency%block(nb)%q_dt

!--- diagnostic tracers are being updated in-place
!--- tracer fields must be returned to the Atm structure
     q(ibs:ibe,jbs:jbe,:,nt_prog+1:ncnst) = Physics_tendency%block(nb)%qdiag

   enddo

   call update_fv_phys( dt_atmos, nt_prog, &
         .true., do_atmos_nudge, Time_next)

!---- diagnostics for FV dynamics -----
   call timing_on('FV_DIAG')
   fv_time = Time + Time_step_atmos
   call get_time (fv_time, seconds,  days)

   call fv_diag(fv_time, nlon, mlat, nlev, beglat, endlat, &
        ncnst, zvir, dt_atmos, .false.)
   call timing_off('FV_DIAG')

 end subroutine atmosphere_state_update


 subroutine atmos_physics_driver_inputs (Physics, Atm_block, Physics_tendency)
#include "fv_arrays.h"
   type (physics_type),  intent(inout) :: Physics
   type (block_control_type), intent(in) :: Atm_block
   type (physics_tendency_type), intent(inout), optional :: Physics_tendency
!--- local variabls
   integer :: nb, ibs, ibe, jbs, jbe
#include "fv_point.inc"

!---------------------------------------------------------------------
! use most up to date atmospheric properties when running serially
!---------------------------------------------------------------------

!$OMP parallel do default(shared) private(nb, ibs, ibe, jbs, jbe)
   do nb = 1, Atm_block%nblks
     ibs = Atm_block%ibs(nb)
     ibe = Atm_block%ibe(nb)
     jbs = Atm_block%jbs(nb)
     jbe = Atm_block%jbe(nb)

     Physics%block(nb)%phis = phis(ibs:ibe,jbs:jbe)
     Physics%block(nb)%u    = ua(ibs:ibe,jbs:jbe,:)
     Physics%block(nb)%v    = va(ibs:ibe,jbs:jbe,:)
     Physics%block(nb)%t    = pt(ibs:ibe,jbs:jbe,:)
     Physics%block(nb)%q    = q(ibs:ibe,jbs:jbe,:,1:nt_prog)
     Physics%block(nb)%omega= omga(ibs:ibe,jbs:jbe,:)
     Physics%block(nb)%pe   = pe(ibs:ibe,:,jbs:jbe)
     Physics%block(nb)%peln = peln(ibs:ibe,:,jbs:jbe)
     Physics%block(nb)%delp = delp(ibs:ibe,jbs:jbe,:)
     if (.not.Physics%control%phys_hydrostatic) &
        call mpp_error(FATAL,'atmosphere: the non-hydrostatic option is not supported (phys_hydrostatic=.false.)')
     if (_ALLOCATED(Physics%block(nb)%tmp_4d)) &
        Physics%block(nb)%tmp_4d = q(ibs:ibe,jbs:jbe,:,nt_prog+1:ncnst)

     call fv_compute_p_z (Atm_block%npz, Physics%block(nb)%phis, Physics%block(nb)%pe, &
                          Physics%block(nb)%peln, Physics%block(nb)%delp, Physics%block(nb)%delz, &
                          Physics%block(nb)%t, Physics%block(nb)%q(:,:,:,Physics%control%sphum), &
                          Physics%block(nb)%p_full, Physics%block(nb)%p_half, &
                          Physics%block(nb)%z_full, Physics%block(nb)%z_half, &
                          Physics%control%phys_hydrostatic)

     if (PRESENT(Physics_tendency)) then
!--- copy the dynamics tendencies into the physics tendencies
!--- if one wants to run physics concurrent with dynamics,
!--- these values would be zeroed out and accumulated
!--- in the atmosphere_state_update

       Physics_tendency%block(nb)%u_dt = 0.
       Physics_tendency%block(nb)%v_dt = 0.
       Physics_tendency%block(nb)%t_dt = 0.
       Physics_tendency%block(nb)%q_dt = 0.
       Physics_tendency%block(nb)%qdiag = q(ibs:ibe,jbs:jbe,:,nt_prog+1:ncnst)
     endif
   enddo

 end subroutine atmos_physics_driver_inputs


 subroutine atmos_radiation_driver_inputs (Time, Radiation, Atm_block )
#include "fv_arrays.h"
   type (time_type),      intent(in)    :: Time
   type (radiation_type), intent(inout) :: Radiation
   type (block_control_type), intent(in) :: Atm_block
!--- local variables
   integer :: nb, ibs, ibe, jbs, jbe
#include "fv_point.inc"

!---------------------------------------------------------------------
! use most up to date atmospheric properties when running serially
!---------------------------------------------------------------------
!$OMP parallel do default(shared) private(nb, ibs, ibe, jbs, jbe)
   do nb = 1,Atm_block%nblks
     ibs = Atm_block%ibs(nb)
     ibe = Atm_block%ibe(nb)
     jbs = Atm_block%jbs(nb)
     jbe = Atm_block%jbe(nb)

     Radiation%block(nb)%phis = phis(ibs:ibe,jbs:jbe)
     Radiation%block(nb)%t    = pt(ibs:ibe,jbs:jbe,:)
     Radiation%block(nb)%q    = q(ibs:ibe,jbs:jbe,:,1:nt_prog)
     Radiation%block(nb)%pe   = pe(ibs:ibe,:,jbs:jbe)
     Radiation%block(nb)%peln = peln(ibs:ibe,:,jbs:jbe)
     Radiation%block(nb)%delp = delp(ibs:ibe,jbs:jbe,:)
     if (.not.Radiation%control%phys_hydrostatic) &
        call mpp_error(FATAL,'atmosphere: the non-hydrostatic option is not supported (phys_hydrostatic=.false.)')

     call fv_compute_p_z (Atm_block%npz, Radiation%block(nb)%phis, Radiation%block(nb)%pe, &
                          Radiation%block(nb)%peln, Radiation%block(nb)%delp, Radiation%block(nb)%delz, &
                          Radiation%block(nb)%t, Radiation%block(nb)%q(:,:,:,Radiation%control%sphum), &
                          Radiation%block(nb)%p_full, Radiation%block(nb)%p_half, &
                          Radiation%block(nb)%z_full, Radiation%block(nb)%z_half, &
                          Radiation%control%phys_hydrostatic)
   enddo

!----------------------------------------------------------------------
! obtain pressure-weighted global mean co2 dry volume mixing ratio for
! use by radiation package.
!----------------------------------------------------------------------
! compute_g_avg must be called here because it contains
! mpp_sums that cannot be called during the concurrent radiation
! phase due to the way in which MPI interacts with nested OpenMP
!----------------------------------------------------------------------
!--- esm2_bugs is a logical switch (default=.FALSE.) read in via the nml
!--- in order to allow early esm2 simulations to reproduce a bug
!--- inside of compute_g_avg
    call compute_g_avg(Time, 'co2', Radiation, Atm_block, esm2_bugs)

 end subroutine atmos_radiation_driver_inputs


 subroutine fv_compute_p_z (npz, phis, pe, peln, delp, delz, pt, q_sph, p_full, p_half, z_full, z_half, hydrostatic)
    integer, intent(in)  :: npz
    real, dimension(:,:),   intent(in)  :: phis
    real, dimension(:,:,:), intent(in)  :: pe, peln, delp, delz, pt, q_sph
    real, dimension(:,:,:), intent(out) :: p_full, p_half, z_full, z_half
    logical, intent(in)  :: hydrostatic
!--- local variables
    integer i,j,k
    real    tvm
    real    :: zvir, rrg, ginv

    zvir = rvgas/rdgas - 1.
    ginv = 1./ grav
    rrg  = rdgas / grav

!----------------------------------------------------
! Compute pressure and height at full and half levels
!----------------------------------------------------
    z_half(:,:,npz+1) = phis(:,:) * ginv

    do k=1,npz+1
      do j=1,size(phis,2)
         do i=1,size(phis,1)
           p_half(i,j,k) = pe(i,k,j)
        enddo
      enddo
    enddo

!--------- Hydrostatic option ----------------------------------------------
    if (hydrostatic ) then
      do k=npz,1,-1
        do j=1,size(phis,2)
          do i=1,size(phis,1)
            tvm = rrg*pt(i,j,k)*(1.+zvir*q_sph(i,j,k))
            p_full(i,j,k) = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
            z_full(i,j,k) = z_half(i,j,k+1) + tvm*(1.-p_half(i,j,k)/p_full(i,j,k))
            z_half(i,j,k) = z_half(i,j,k+1) + tvm*(peln(i,k+1,j)-peln(i,k,j))
          enddo
        enddo
      enddo
    else
!--------- Non-Hydrostatic option ------------------------------------------
      do k=npz,1,-1
        do j=1,size(phis,2)
          do i=1,size(phis,1)
            p_full(i,j,k) = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
            z_half(i,j,k) = z_half(i,j,k+1) - delz(i,j,k)
            z_full(i,j,k) = 0.5*(z_half(i,j,k) + z_half(i,j,k+1))
          enddo
        enddo
      enddo
    endif

 end subroutine fv_compute_p_z

 subroutine reset_atmos_tracers (Physics, Physics_tendency, Atm_block)
#include "fv_arrays.h"
   type (physics_type), intent(in) :: Physics
   type (physics_tendency_type), intent(in) :: Physics_tendency
   type (block_control_type), intent(in) :: Atm_block
!--- local variables
   integer :: nb, ibs, ibe, jbs, jbe
#include "fv_point.inc"

!--- After initialization by the physics, tracer fields must be
!--- returned to the Atm structure.  This is because tracer driver
!--- can reset the initial values
!$OMP parallel do default(shared) private(nb, ibs, ibe, jbs, jbe)
    do nb = 1, Atm_block%nblks
      ibs = Atm_block%ibs(nb)
      ibe = Atm_block%ibe(nb)
      jbs = Atm_block%jbs(nb)
      jbe = Atm_block%jbe(nb)

      q(ibs:ibe,jbs:jbe,:,1:nt_prog)       = Physics%block(nb)%q
      q(ibs:ibe,jbs:jbe,:,nt_prog+1:ncnst) = Physics_tendency%block(nb)%qdiag
    enddo

 end subroutine reset_atmos_tracers

end module atmosphere_mod
