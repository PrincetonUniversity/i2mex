!
! ------------------- tree_close ------------------
! close the open mdsplus tree
!
subroutine tree_close(ilun, ishot, ztree, ier)
  implicit none
 
  integer,       intent(in)  :: ilun   ! logical unit for messages
  integer,       intent(in)  :: ishot  ! shot number
  character*(*), intent(in)  :: ztree  ! tree name
  integer,       intent(out) :: ier    ! error flag
 
  integer,external :: mds_close, mds_value, idescr_long
  integer :: istat,stat,retl
 
  ier=0
  !istat = mds_close(ztree, ishot)
  stat=mds_value('tcl("close")',idescr_long(istat),retl)
  if (stat.ne.1) then
     write(ilun,'(/a)') '?trdat_tomds.tree_close: error closing tree '//ztree
     ! call mdserr(ilun, 'tree close status:  ',istat)
     ier=1
  end if
end subroutine tree_close
 
