subroutine trx_label(zlbl)
 
  use trx_module
  implicit NONE
 
  character*(*) zlbl
 
!
!  return label for current TRANSP run
!  if no run is connected, a blank label is returned
!
  if(nsurf.eq.0) then
     zlbl=' '
  else
     zlbl=run_label
  endif
!
  return
end subroutine trx_label
