 
!
! ------------------- wr_longarr -----------------
!
! put a long integer array into the tree
!
subroutine wr_longarr(ilun, mds_path, name, value, nvalue, ier)
 
  implicit none
 
  integer,       intent(in)    :: ilun         ! logical unit for messages
  integer,       intent(in)    :: nvalue         ! size of value
  character*(*), intent(in)    :: mds_path     ! path to current node
  character*(*), intent(in)    :: name         ! new signal node to create
  integer, dimension(:), intent(in)    :: value(nvalue)        ! string to place in tree
  integer,       intent(out)   :: ier          ! error flag
 
  logical,external :: lmds_errstat
  integer,external :: MdsPut, idescr_longarr, mds_add_node
  integer          :: dsc, stat, nid, n, idims(1)
  character*150    :: expr
 
  ier=0
  expr = trim(mds_path)//':'//trim(name)
 
  print *, '%trdatbuf_tomds: putting ' // trim(expr)
 
  nid = mds_add_node(trim(expr),'numeric',stat)
  if (nid .le. 0) then
     write(ilun,'(/a)') '?trdatbuf_tomds.wr_longarr: error creating node '//expr
     ier=1
     return
  end if
 
  n=1
  idims(1)=nvalue
  dsc = idescr_longarr(value, idims, n)
  stat=MdsPut(trim(expr)//char(0),'$'//char(0),dsc,0)
  if (lmds_errstat(ilun, stat)) then
     write(ilun,'(/a)') '?trdatbuf_tomds.wr_longarr: error putting '//name
     ier=1
  end if
 
end subroutine wr_longarr
 
