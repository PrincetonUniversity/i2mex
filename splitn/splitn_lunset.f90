subroutine splitn_lunset(ilun)

  use splitn_module
  implicit NONE

  ! set FORTRAN LUN for splitn i/o operations

  integer, intent(in) :: ilun

  lun = ilun

end subroutine splitn_lunset

subroutine splitn_lunget(jlun)

  use splitn_module
  implicit NONE

  ! get FORTRAN LUN for splitn i/o operations

  integer, intent(out) :: jlun

  jlun = lun

end subroutine splitn_lunget
