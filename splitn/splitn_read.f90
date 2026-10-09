subroutine splitn_read(fname,ios)

  !  READ a TRANSP namelist into the splitn_module ...

  use splitn_module
  implicit NONE

  character*(*), intent(in) :: fname     ! filename
  integer, intent(out) :: ios            ! completion status code, 0=OK

  !--------------------------------------

  call read_nl(fname,ios)

end subroutine splitn_read
