!
! saving the server can cause problems if mdsplus accesses a different
! server through other software
!
subroutine mds_clear_save
  use cplotr_mod
  implicit none
  integer is, MdsSetSocket  

  mds_save_server = ' '
  is = MdsSetSocket(-1)
end subroutine mds_clear_save
