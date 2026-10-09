 
!
! ------------------- wr_cstring -----------------
!
! put a character string into the tree
!
subroutine wr_cstring(ilun, mds_path, name, value, ier)
 
  implicit none
 
  integer,       intent(in)    :: ilun         ! logical unit for messages
  character*(*), intent(in)    :: mds_path     ! path to current node
  character*(*), intent(in)    :: name         ! new signal node to create
  character*(*), intent(in)    :: value        ! string to place in tree
  integer,       intent(out)   :: ier          ! error flag
 
  logical,external :: lmds_errstat
  integer,external :: MdsPut, idescr_cstring, mds_add_node
  integer          :: dsc, stat, nid
  character*150    :: expr
 
  ier=0
  expr = trim(mds_path)//':'//trim(name)
 
  print *, '%trdatbuf_tomds: putting ' // trim(expr)
 
  nid = mds_add_node(trim(expr),'text',stat)
  if (nid .le. 0) then
     write(ilun,'(/a)') '?trdatbuf_tomds.wr_cstring: error creating node '//expr
     ier=1
     return
  end if
 
  dsc = idescr_cstring(value)
  stat=MdsPut(trim(expr)//char(0),'$'//char(0),dsc,0)
  if (lmds_errstat(ilun, stat)) then
     write(ilun,'(/a)') '?trdatbuf_tomds.wr_cstring: error putting '//name
     ier=1
  end if
 
end subroutine wr_cstring
 
