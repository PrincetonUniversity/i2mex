logical function tdb_steerxcep(id,ilev)
  !------------------------------------------------
  ! generally, trdatbuf identifiers beginning with L (addresses)
  ! are steerable; those beginning with other letters (usually N) are
  ! not, with the exception of identifiers defined in the DATA statement
  ! here...
  !------------------------------------------------
  character*(*), intent(in) :: id
  integer,intent(out) :: ilev    ! see explanation...

  ! on output, if tdb_steerxcep = .TRUE.:
  !   ilev=0 means the item value can change but only if both instances are
  !          non-zero
  !   ilev=1 means the item value can change with no restriction
  !
  !------------------------------------------------
  integer :: i
  character*16 :: zid
  integer, parameter :: nids = 15
  character*16 :: zid_xcep(nids)

  data zid_xcep/'NTIME1','NTIME2', &
       '~NTSAW','~NSAWFLAG','~LPELDA','~NPELDA', &
       'NTIMNB','NTIMRF','NTIMRFF','NTIMEC','NTIMLH', &
       'NDMMX','NDFS','NTPSI','NDPSI'/
  !------------------------------------------------

  zid=id
  call uupper(zid)

  ilev = -1

  if(zid(1:1).eq.'L') then
     tdb_steerxcep = .TRUE.
     ilev = 0
     do i=1,nids
        if('~'//trim(zid).eq.trim(zid_xcep(i))) ilev = 1
     enddo

  else
     tdb_steerxcep = .FALSE.
     do i=1,nids
        if(trim(zid).eq.trim(zid_xcep(i))) then
           tdb_steerxcep=.TRUE.
           ilev = 0
        else if('~'//trim(zid).eq.trim(zid_xcep(i))) then
           tdb_steerxcep=.TRUE.
           ilev = 1
        endif
     enddo
  endif

end function tdb_steerxcep
