subroutine trx_bxtr_method_set(imethod)

  use trx_bxtr_options
  implicit NONE

  integer, intent(in) :: imethod
  !  imethod = 1 -- ad hoc extrapolator (recommended & default)
  !  imethod = 2 -- del.B = 0 extrapolator

  extrap_method = max(1,min(2,imethod))

end subroutine trx_bxtr_method_set

subroutine trx_bxtr_method_get(imethod)

  use trx_bxtr_options
  implicit NONE

  integer, intent(out) :: imethod
  !  imethod = 1 -- ad hoc extrapolator (recommended & default)
  !  imethod = 2 -- del.B = 0 extrapolator

  imethod = extrap_method

end subroutine trx_bxtr_method_get

subroutine trx_bxtr_smooth_set(zsmooth)

  use trx_bxtr_options
  implicit NONE

  real*8, intent(in) :: zsmooth ! edge smoothing parameter (meters)
  real*8, parameter :: ZERO = 0.0d0

  edge_smooth = max(ZERO,zsmooth)

end subroutine trx_bxtr_smooth_set

subroutine trx_bxtr_smooth_get(zsmooth)

  use trx_bxtr_options
  implicit NONE

  real*8, intent(out) :: zsmooth ! edge smoothing parameter (meters)

  zsmooth = edge_smooth

end subroutine trx_bxtr_smooth_get

