subroutine trx_rzbc_fdiff

  use trx_module
  implicit NONE

  isign_rzbc=-1   ! use finite difference BC on R,Z equilibrium profiles

end subroutine trx_rzbc_fdiff

subroutine trx_rzbc_nknot

  use trx_module
  implicit NONE

  isign_rzbc=1    ! use not a knot BC on R,Z equilibrium profiles

end subroutine trx_rzbc_nknot
