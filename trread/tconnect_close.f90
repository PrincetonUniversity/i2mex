subroutine tconnect_close(ierr)
 
  use tconnect_mod
 
  implicit NONE
 
  integer ierr
 
  ! -----------------------------------------
 
  if(tc_mds_open.eq.1) then
     call mdscls(tc_idrun,ierr)
     tc_mds_open=0
  endif
 
  return
 
end subroutine tconnect_close
