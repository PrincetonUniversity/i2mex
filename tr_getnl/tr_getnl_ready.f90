subroutine tr_getnl_ready(iflag)

  !  return iflag = .TRUE. if the tr_getnl module contains a namelist

  use tr_getnl
  implicit NONE

  logical, intent(out) :: iflag



  iflag = (mltext_nlines.gt.0)

end subroutine tr_getnl_ready
