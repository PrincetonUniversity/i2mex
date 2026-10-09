!
! ----------------- tclwrite ------------------
! send WRITE command to mdsplus
!
subroutine tclwrite(ilun, ier)
  implicit none
 
  integer,       intent(in)  :: ilun   ! logical unit for messages
  integer,       intent(out) :: ier    ! error flag
 
  integer           :: stat,istat,retl
  integer, external :: mds_value, idescr_long
 
  ier=0
  write(ilun,*) '%tclwrite: writing to tree'
  stat=mds_value('tcl("write")',idescr_long(istat),retl)
  if(stat.ne.1) then
     write(ilun,'(/a)') '?trdat_tomds.tclwrite: error in tcl write'
     call mdserr(ilun, '+tclwrite status:  ',istat)
     ier=1
  endif
end subroutine tclwrite
 
!
! ----------------- tcldelete ------------------
! send DELETE command to mdsplus
!
subroutine tcldelete(ilun, ier)
  implicit none
 
  integer,       intent(in)  :: ilun   ! logical unit for messages
  integer,       intent(out) :: ier    ! error flag
 
  integer           :: stat,istat,retl
  integer, external :: mds_value, idescr_long
 
  ier=0
  write(ilun,*) '%tcldelete: deleting .TRDATA node'
  stat=mds_value('tcl("delete node .TRDATA /noconfirm")',idescr_long(istat),retl)
  if(stat.ne.1) then
     write(ilun,'(/a)') '?trdat_tomds.tcldelete: error in tcl delete'
     call mdserr(ilun, '+tclwrite status:  ',istat)
     ier=1
  endif
end subroutine tcldelete
