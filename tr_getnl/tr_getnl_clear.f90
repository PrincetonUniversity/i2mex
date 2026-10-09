subroutine tr_getnl_clear
!
!  clear the namelist module
!
  use tr_getnl
  implicit NONE
!
!  OK... init namelist status (read not yet tried)
!
  nltext_status=-1
  nltext_nlines=0
  if(allocated(nltext)) deallocate(nltext)
  if(allocated(nltext_lens)) deallocate(nltext_lens)
!
  return
!
end subroutine tr_getnl_clear
