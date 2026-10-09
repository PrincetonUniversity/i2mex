 
!
! ------------------- wr_cstringarr -----------------
!
! put a character string array into the tree
!
subroutine wr_cstringarr(ilun, mds_path, name, value, nvalue, ier)
 
  implicit none
 
  integer,       intent(in)    :: ilun         ! logical unit for messages
  integer,       intent(in)    :: nvalue         ! size of value
  character*(*), intent(in)    :: mds_path     ! path to current node
  character*(*), intent(in)    :: name         ! new signal node to create
  character*(*), intent(in)    :: value(nvalue)        ! string to place in tree
  integer,       intent(out)   :: ier          ! error flag
 
  logical,external :: lmds_errstat
  integer,external :: mdsPut, idescr_cstringarr, mds_add_node
  integer          :: dsc, stat, nid, n, idims(1)
  character*150    :: expr
 
  ier=0
  expr = trim(mds_path)//':'//trim(name)
 
  print *, '%trdatbuf_tomds: putting ' // trim(expr)
 
  nid = mds_add_node(trim(expr),'text',stat)
  if (nid .le. 0) then
     write(ilun,'(/a)') '?trdatbuf_tomds.wr_cstringarr: error creating node '//expr
     ier=1
     return
  end if
 
  n=1
  idims(1)=nvalue
  dsc = idescr_cstringarr(value, idims, n)
  stat=MdsPut(trim(expr)//char(0),'$'//char(0),dsc,0)
  if (lmds_errstat(ilun, stat)) then
     write(ilun,'(/a)') '?trdatbuf_tomds.wr_cstringarr: error putting '//name
     ier=1
  end if
 
end subroutine wr_cstringarr
 
