subroutine f_upwind0(UPWIND0,zdiffus,zdravi,zveloc,zphip1,zphi)

  use iso_c_binding, only: fp => c_double
  use physconst_mod, only: zero, half, one
  implicit none

  !  evaluate upwind adjustment
  !  on output:  zphip1=0.5, zphi=0.5 => no upwind adjustment
  !  on output:  zphip1=0.0, zphi=1.0 => full upwind adjustment (outflow)
  !  on output:  zphip1=1.0, zphi=0.0 => full upwind adjustment (inflow)
  !    intermediate adjustments possible; zphi+zphip1 = 1.0 always.

  !-------------------------------------------

  real(fp), intent(in) :: UPWIND0   ! upwind parameter (dimensionless)

  !  roughly, this limits fluctuations df/f that can occur due to lack
  !  of upwind differencing; the smaller the value set, the stronger
  !  the upwind adjustment that will be made

  !  the following dimensional quanties can be cgs or mks but choose
  !  one or the other

  real(fp), intent(in) :: zdiffus   ! diffusivity  (cm2/sec) or (m2/sec)  .ge. 0
  real(fp), intent(in) :: zdravi    ! <1/dr> (1/cm) or (1/m) -- surface spacing
  real(fp), intent(in) :: zveloc    ! flow velocity (cm/sec) or (m/sec)

  !  output

  real(fp), intent(out) :: zphip1   ! weight on outward zone
  real(fp), intent(out) :: zphi     ! weight on inward zone

  !---------------------------------------
  ! local:
  real(fp) :: zalph
  !---------------------------------------

  if(zdiffus.eq.zero) then

    !  special case: DIFFUSIVITY = ZERO

    if(zveloc.gt.zero) then
      zphi=one
      zphip1=0.0

    else if(zveloc.lt.zero) then
      zphi=0.0
      zphip1=one

    else
      zphi=half
      zphip1=half
    endif

    return

  endif

  if(zveloc.eq.zero) then

    !  special case: VELOCITY = ZERO

    zphi=half
    zphip1=half

    return

  endif

  !  normal case:  both D and v are non-zero; look at ratio...
  !  normalizing factor UPWIND0 from TRANSP namelist...

  zalph = min(one, UPWIND0*zdiffus*zdravi/abs(zveloc) )

  if(zveloc.gt.zero) then

    zphip1=zalph*half
    zphi=one-zphip1

  else

    zphi=zalph*half
    zphip1=one-zphi

  endif

  return
end subroutine f_upwind0
