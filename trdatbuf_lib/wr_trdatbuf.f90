subroutine wr_trdatbuf(nonlin, d, zserver, ztree, zshot, ztok, iyear, ishot, ier)
  Use trdatbuf_module
  Use trdatbuf_intmod

  implicit none


  integer :: ishot       ! mds+ shot
  integer :: ier         ! error flag
  integer :: isocket     ! from conopn
  integer :: iyear       ! two digit year
  integer :: nonlin      ! lun for error messages
  integer :: ilen        !  size of datbuf
  integer :: iparams(3)     !  sizes of MDSplus arrays

  character*80  :: zshot    ! runid
  character*120 :: zserver  ! mds+ server
  character*80  :: ztree    ! mds+ tree
  character*20  :: ztok     ! tokamak
  character*10  :: zyear    ! year of shot
  character*10  :: zpulse   ! pulse string
  character*80  :: zerr     ! error string
  character*40  :: zchild   ! temp
  character*1   :: bcks

  integer :: stat

  Type (trdatbuf) :: d

  call mds_conopn(nonlin, zserver, ztree, zshot, ztok, iyear, ishot, isocket, 4, ier)
  if (ier/=0) then
     write(nonlin,'(/a)')     '?trdatbuf_tomds: error opening the tree'
     goto 999
  else
!     write(*,*) ' mds_conopn successful'
  end if

! 10/11/05 CLF: don't delete
!   call tcldelete(nonlin,ier)  ! delete .TRBUFDATA node

  !
  ! ------------- start processing -----------------
  !
  !
  ! ------------- build .TRBUFDATA child ---------------
  ! Create .TRBUFDATA child
  !
  zchild = '.TRBUFDATA'
  call trnode_child(nonlin, zchild, ier)
  if (ier/=0) goto 990
!  write(*,*) ' trnode_child2 successful'

!
! write DATBUF real data

  ILEN = d%LFREE - 1
  call wr_doublearr(nonlin, zchild, 'DATBUF', d%datbuf, ILEN, ier)
  if (ier/=0) goto 990
!  write(*,*) ' wr_doublearr successful '

  call intmod_init
  call fill_intbuf(d)
!
! write IVARSIZE integer data
  call wr_longarr(nonlin, zchild, 'IVARSIZE', ivarsize, nivar, ier)
  if (ier/=0) goto 990
!  write(*,*) ' wr_longarr successful '

!
! write IVARNAME character data
  call wr_cstringarr(nonlin, zchild, 'IVARNAME', ivarname, nivar, ier)
  if (ier/=0) goto 990
!  write(*,*) ' wr_cstringarr successful '
!
! write INTBUF integer data
  ilen2=size(intbuf)
  call wr_longarr(nonlin, zchild, 'INTBUF', intbuf, ilen2, ier)
  if (ier/=0) goto 990
!  write(*,*) ' wr_longarr successful '

  iparams = (/nivar, ilen, ilen2/)
!  write(*,*) "iparams=",iparams
!
! write PARAMS integer data
  call wr_longarr(nonlin, zchild, 'PARAMS', iparams, 3, ier)
  if (ier/=0) goto 990
!  write(*,*) ' wr_longarr successful '


  call tclwrite(nonlin, ier)                 ; if (ier/=0) goto 990
  call tree_close(nonlin, ishot, ztree, ier) ; if (ier/=0) goto 990
!  write(*,*) "tree_close2 OK"

990 continue
  call MdsCacheDisconnect
  stat=0
999 continue


  return
end subroutine wr_trdatbuf
