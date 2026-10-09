subroutine splitn_addwarn(item0,zlist,nlist,iwid,klist,maxlist)

  ! build alphabetic list; discard item without warning if it is a duplicate
  ! or if the maximum list size has been reached.
  !
  ! these are expected to be very short lists so we don't bother with n*log(n)
  ! sorting methods.

  use splitn_module, only: kupdate
  implicit NONE

  character*(*), intent(in) :: item0 ! item to add
  integer, intent(inout) :: nlist    ! current list size (prior to add)
  integer, intent(in) :: iwid        ! width of list elements
  integer, intent(in) :: maxlist     ! maximum list size
  character*(iwid), intent(inout) :: zlist(maxlist)  ! list being added to
  integer, intent(inout) :: klist(maxlist)  ! update block index
  
  !-------------------------------------------

  character*(iwid) item
  integer i,ipos

  !-------------------------------------------

  item=item0

  if(nlist.eq.maxlist) return         ! ignore max length exceeded

  ipos=nlist+1
  do i=1,nlist
     if((item.eq.zlist(i)).and.(kupdate.eq.klist(i))) then
        return  ! ignore duplicates
     endif
     if(item.le.zlist(nlist)) ipos=i
  enddo

  !  insert in list, maintaining alphabetic order

  do i=nlist,ipos,-1
     zlist(i+1)=zlist(i)
     klist(i+1)=klist(i)
  enddo

  zlist(ipos)=item
  klist(ipos)=kupdate

  nlist=nlist+1

end subroutine splitn_addwarn
