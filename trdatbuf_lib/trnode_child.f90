!
! ------------------- trnode_child ---------------------
! create a child node
!
subroutine trnode_child(ilun, mds_path, ier)
  implicit none
 
  integer,       intent(in)  :: ilun   ! logical unit for messages
  character*(*), intent(in)  :: mds_path   ! path to child node
  integer,       intent(out) :: ier        ! error flag
 
  integer,external :: mds_add_node
  integer stat,nid
 
  ier=0
 
  print *, '%trnode_child: putting ' // trim(mds_path)
 
  nid=mds_add_node(trim(mds_path),'structure',stat)
  if (nid.le.0) then
     if (stat .ne. 265388168) then  ! TreeALREADY_THERE
        write(ilun,'(/a)') '?trnode_child: error creating node '//mds_path
        call mdserr(ilun, trim(mds_path)//' node create status:  ',stat)
        ier=1
     end if
  end if
end subroutine trnode_child
 
