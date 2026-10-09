subroutine trx_init
 
! initialize/clear trx -- transp data access & interpolation
! ...use when "connecting" to a new run
 
  use trx_module
  implicit NONE
 
  INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)

  integer :: init = 0
  save init
 
! -------------------------
 
  d_data_avail = .FALSE.

  !  d_init_flag is initialized once only in the module
  !  (should not be reset)

  if(init.eq.0) then
     call initcpl
     init=1
  endif
 
  nsctime=0
  nprtime=0
  nxmax=0
  nmax=0
  nsurf=0
  nrmajm=0
  nrmjsym=0
 
  ksym=0
  nmoms=0

  itype_x=0
  itype_xb=0
  itype_Rmajm=0
  itype_Rmjsym=0

  id_rho=0
  id_rhozc=0
  id_chi=0
  id_R=0
  id_Z=0
  id_Rmajm=0
  id_Rmjsym=0

  run_label = ' '
  tdev = ' '

  tmin=0.0_R8
  tmax=0.0_R8
 
  itset=0
  time0=0.0_R8
  delta_t=0.0_R8
 
  call trx_bc_init
 
  units_status= -1
  xaxis_id = -1
  prof_iord = -1
 
  call tr_getnl_clear

!------------------------------

  mds_arch_flag = .FALSE.
  file_path = ' '

!------------------------------

  nRfree=0
  nZfree=0
  ifound_Psi0=0
  if(allocated(Rgrid_free)) deallocate(Rgrid_free)
  if(allocated(Zgrid_free)) deallocate(Zgrid_free)
  if(allocated(PsiRZ_free)) deallocate(PsiRZ_free)
  raxis_mhd=0.0_R8
  zaxis_mhd=0.0_R8
  Psi0_mhd=0.0_R8

  return
 
end subroutine trx_init
 
subroutine trx_bc_init
 
! (re)initialize standard boundary conditions & delete over-rides
 
  use trx_module
  implicit NONE
 
! --------------------
!
!  defaults:  f' -> 0 on axis
!  unconstrained at edge
!
  INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
!
  ibc_axis = 1
  zbc_axis = 0.0E0_R8
  ibc_edge = 0
  zbc_edge = 0.0_R8    ! (unused unless ibc_edge is changed)
!
!  defaults:  no over-rides
!
  ival_axis = 0
  zval_axis = 0.0_R8   ! (unused unless an over-ride is set)
  ival_edge = 0
  zval_edge = 0.0_R8   ! (unused unless an over-ride is set)
!
!  defaults:  no sign changing extrapolation
!
  isign_check = 1
  zsign_lim =0.25_R8   ! extrapolation not below 1/4 * sign(Edge)*mod(Edge)
  return
!
end subroutine trx_bc_init
