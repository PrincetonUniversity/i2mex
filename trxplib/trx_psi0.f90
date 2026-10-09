subroutine trx_psi0(ifind,psi0)

  !  return Psi0 data (poloidal flux diff, mag. axis to machine axis)
  !  if available

  use trx_module
  implicit NONE

  !  result is based on most recently read TRANSP data

  integer, intent(out) :: ifind   ! =1: data available; =0: not available
  real*8,  intent(out) :: Psi0    ! Psi0 data (if ifind=1) or ZERO (if ifind=0)

  !-----------------------

  if(ifound_psi0.eq.1) then
     ifind=1
     Psi0 = Psi0_mhd
  else
     ifind=0
     Psi0 = 0.0d0
  endif

end subroutine trx_psi0
