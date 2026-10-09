subroutine tdb_find_next_sawtime(d,zprev_time,ztime)

  ! get the next sawtooth time -- a time greater than zprev_time but
  ! also between TINIT and FTIME.  If no such time exists, a very large
  ! value us returned for ztime; all times in SECONDS

  use trdatbuf_obj
  implicit NONE

  type (trdatbuf) :: d
  real*8, intent(in) :: zprev_time   ! preceding time (input)
  real*8, intent(inout) :: ztime     ! time found or a huge time (output)

  !----------
  integer :: it
  real*8 :: ztnext
  real*8, parameter :: zlarge = 1.0D34
  !----------

  if(d%ltsaw.eq.0) then
     ztime = zlarge      ! in this case there are no sawteeth, ever
  else if (d%datbuf(d%ltsaw+d%ntsaw-1).le.zprev_time) then
     ztime = zlarge      ! no more sawteeth
  else
     ! find the next sawtooth.  This is not an efficient search.  It can be
     ! improved should that ever be necessary.
     do it = 1,d%ntsaw
        ztnext = d%datbuf(d%ltsaw+it-1)
        if(ztnext.gt.zprev_time) then
           ztime = ztnext
           exit
        endif
     enddo
  endif

end subroutine tdb_find_next_sawtime

subroutine tdb_num_sawtimes(d,ntsaw)
  use trdatbuf_obj
  implicit NONE

  type (trdatbuf) :: d
  integer, intent(out) :: ntsaw  ! number of sawtooth times in trdat buffer

  ntsaw = d%ntsaw
end subroutine tdb_num_sawtimes

subroutine tdb_sawtimes(d,ztimes,inum)
  use trdatbuf_obj
  implicit NONE

  !  return sawtooth times; print warning if inum does not match the
  !  actual number.  If inum is too small, the first inum sawtooth times
  !  are returned; if inum is too large, the extra array elements are set
  !  to a large number.

  type (trdatbuf) :: d
  integer, intent(in) :: inum
  real*8, intent(out) :: ztimes(inum)

  integer :: inuma,lunmsg_tdb
  real*8, parameter :: zlarge = 1.0d34

  if(inum.ne.d%ntsaw) then
     write(lunmsg_tdb(0),*) ' %tdb_sawtimes warning: argument inum=',inum
     write(lunmsg_tdb(0),*) '  this does not match the number of sawtooth'
     write(lunmsg_tdb(0),*) '  event times in the trdat buffer: ',d%ntsaw
     if(inum.gt.d%ntsaw) then
        ztimes(d%ntsaw+1:inum)=zlarge
     endif
  endif

  inuma=min(inum,d%ntsaw)
  ztimes(1:inuma)=d%datbuf(d%ltsaw:d%ltsaw+inuma-1)
end subroutine tdb_sawtimes

subroutine tdb_post_sawtime(d,ztime_pre,ztime_post,iwarn)

  !  given event start (SAWTOOTH) time return event end time

  use trdatbuf_obj
  implicit NONE

  type (trdatbuf) :: d
  real*8, intent(in) :: ztime_pre   ! time at start of event
  real*8, intent(out) :: ztime_post ! time at end of event
  integer, intent(out) :: iwarn     ! 0=OK; 1=no sawtooth events
                                    ! 2= sawtooth events exist but ztime_pre
                                    !    does not match any of them...

  !  if the warning flag is set to a non-zero value, ztime_post = ztime_pre
  !  is returned.

  !---------------------
  real*8 :: zdt_avg, zdt_min, zdt_tol, zdt_fnd
  integer :: ilt,int,ilas,it
  !---------------------

  if(d%nsawflag.ne.1) then
     iwarn=1
     ztime_post=ztime_pre
     return
  endif

  ilt=d%ltime2
  int=d%ntime2
  ilas=ilt+int-1

  zdt_avg=(d%datbuf(ilas) - d%datbuf(ilt))/int
  zdt_min=minval(d%datbuf(ilt+1:ilas)-d%datbuf(ilt:ilas-1))
  !  zdt_min should be the sawtooth event duration
  zdt_tol=zdt_min/5

  iwarn=-1
  do it=1,int-1
     if(abs(ztime_pre-d%datbuf(ilt+it-1)).le.zdt_tol) then
        zdt_fnd=d%datbuf(ilt+it)-d%datbuf(ilt+it-1)
        if(zdt_fnd.gt.zdt_min+zdt_tol) then
           iwarn=2
           ztime_post=ztime_pre
        else
           iwarn=0
           ztime_post=d%datbuf(ilt+it)  ! the next time...
        endif
        exit
     endif
  enddo

  if(iwarn.eq.-1) then
     iwarn=2
     ztime_post=ztime_pre
  endif

end subroutine tdb_post_sawtime

subroutine tdb_post_peltime(d,ztime_pre,ztime_post,iwarn)

  !  given event start (PELLET) time return event end time

  use trdatbuf_obj
  implicit NONE

  type (trdatbuf) :: d
  real*8, intent(in) :: ztime_pre   ! time at start of event
  real*8, intent(out) :: ztime_post ! time at end of event
  integer, intent(out) :: iwarn     ! 0=OK; 
                                    ! 2= no such event time found...

  !  if the warning flag is set to a non-zero value, ztime_post = ztime_pre
  !  is returned.

  !---------------------
  real*8 :: zdt_avg, zdt_min, zdt_tol, zdt_fnd
  integer :: ilt,int,ilas,it
  !---------------------

  ilt=d%ltime2
  int=d%ntime2
  ilas=ilt+int-1

  zdt_avg=(d%datbuf(ilas) - d%datbuf(ilt))/int
  zdt_min=minval(d%datbuf(ilt+1:ilas)-d%datbuf(ilt:ilas-1))
  !  zdt_min should be the pellet event duration
  zdt_tol=zdt_min/5

  iwarn=-1
  do it=1,int-1
     if(abs(ztime_pre-d%datbuf(ilt+it-1)).le.zdt_tol) then
        zdt_fnd=d%datbuf(ilt+it)-d%datbuf(ilt+it-1)
        if(zdt_fnd.gt.zdt_min+zdt_tol) then
           iwarn=2
           ztime_post=ztime_pre
        else
           iwarn=0
           ztime_post=d%datbuf(ilt+it)  ! the next time...
        endif
        exit
     endif
  enddo

  if(iwarn.eq.-1) then
     iwarn=2
     ztime_post=ztime_pre
  endif

end subroutine tdb_post_peltime
