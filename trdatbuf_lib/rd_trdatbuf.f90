subroutine rd_trdatbuf(nonlin, d, zserver, ztree, zshot, ztok, iyear, ishot,  ier)

  !  1.  connect to server & open tree
  !  2.  read trdatbuf data
  !  3.  disconnect

  Use trdatbuf_module

  implicit none

  integer :: ishot       ! mds+ shot
  integer :: ier         ! error flag
  integer :: isocket     ! from conopn
  integer :: iyear       ! two digit year
  integer :: nonlin      ! lun for error messages

  character*80  :: zshot    ! runid
  character*120 :: zserver  ! mds+ server
  character*80  :: ztree    ! mds+ tree
  character*20  :: ztok     ! tokamak
  character*40  :: zchild   ! temp

  integer       :: stat

  Type (trdatbuf) :: d

  call mds_conopn(nonlin, zserver, ztree, zshot, ztok, iyear, ishot, isocket, 3, ier)
  if (ier/=0) then
     write(nonlin,'(/a)')     '?trdatbuf_frommds: error opening the tree'
     goto 999     
  else
!     write(*,*) ' mds_conopn successful'
  end if

  call rd_trdatbuf_mds(d,ier)

  if (ier/=0) then
     write(nonlin,'(/a)')     '?trdatbuf_frommds: error reading the tree'
     go to 990
  endif

990 continue

  call MdsCacheDisconnect
  stat=0

999 continue

  return
end subroutine rd_trdatbuf


subroutine rd_trdatbuf_mds(d,ier)

  !  read trdatbuf data from open tree
  !  do not connect, open, or disconnect.

  Use trdatbuf_module
  Use trdatbuf_intmod

  implicit NONE

  Type (trdatbuf) :: d
  integer, intent(out) :: ier   ! status code returned (0=OK)

  !----------------------
  integer :: ilen        !  size of datbuf
  integer :: nivarx,   ni_varx    ! number of integer variables in ivarname & ivarsize
  integer :: iparams(3)     !  sizes of MDSplus arrays
  integer :: isize       ! size of datbuf returned from datbuf_expand

  real*8, dimension(:), allocatable :: datbuf
  !----------------------

  ier = 0

  call rd_longarr('.TRBUFDATA:PARAMS',iparams,3,ier)
  if (ier/=0) goto 990

  nivarx = iparams(1)
  ilen = iparams(2)
  ilen2 = iparams(3)

  if(allocated(intbuf)) then
     if(size(intbuf).ne.ilen2) then
        deallocate(intbuf)  ! not sure if the size will really change...
        allocate(intbuf(ilen2))
     endif
  else
     allocate(intbuf(ilen2))
  endif

  call rd_longarr('.TRBUFDATA:INTBUF',intbuf,ilen2,ier)
  if (ier/=0) goto 990

  allocate(datbuf(ilen))
  call rd_doublearr('.TRBUFDATA:DATBUF',datbuf, ILEN,ier)
  if (ier/=0) goto 990

  ivarname= ' '  ! clear to blank first...
  ni_varx=size(ivarname)
  call rd_cstringarr('.TRBUFDATA:IVARNAME',ivarname, ni_varx,ier)
  if (ier/=0) goto 990

  call rd_longarr('.TRBUFDATA:IVARSIZE',ivarsize,ni_varx,ier)
  if (ier/=0) goto 990

  call trdatbuf_init(d)

  call get_intbuf(d)
  call datbuf_expand(d,ilen,isize)

  d%DATBUF(1:ilen) = datbuf

990 continue

  return
end subroutine rd_trdatbuf_mds
