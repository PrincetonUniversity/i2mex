subroutine trx_uget(ident,ixid,zname,zunits,istat,ierr)
!
  use trx_module
  implicit NONE
!
!  return the name & units of the item associated with xplasma id "ident".
!
  integer, intent(in) :: ident         ! xplasma id of item
  integer, intent(out) :: ixid         ! x-axis id of item (0 if not 1d).
  character(10), intent(out) :: zname  ! name of item
  character(10), intent(out) :: zunits ! physical units
  integer, intent(out) :: istat        ! MKS status (0=OK, 1=warning)
  integer, intent(out) :: ierr         ! completion code, 0=OK
!
!---------------------------------------
! get xplasma name
!
  call eq_get_fname(ident,zname)
  if(zname.eq.' ') then
     ierr=1
  else
     ierr=0
     zunits=prof_units(ident)
     istat=units_status(ident)
     if(istat.lt.0) ierr=1
     ixid=xaxis_id(ident)
     if(ixid.lt.0) ierr=1
  endif
!
  return
end subroutine trx_uget
 
subroutine trx_uget_next(ident)
!
  use trx_module
  implicit NONE
!
!   return id of next item with units
!
  integer, intent(inout) :: ident      ! search start / search result
!
!   on input, ident gives search start point:  start at (ident+1)
!  on output, ident=N means item N is next; ident=0 means no more items.
!
  integer jdent
!
  do while (ident.lt.maxid)
     ident=ident+1
     if(units_status(ident).ge.0) return
  end do
  ident = 0
  return
end subroutine trx_uget_next
