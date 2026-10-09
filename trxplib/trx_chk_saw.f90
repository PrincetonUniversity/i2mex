subroutine trx_chk_saw_r8(ztime,ievent)

  ! check time (call single precision routine)
  ! see trx_chk_saw comments

  implicit NONE

  real*8, intent(inout) :: ztime ! time to (possibly) be adjusted
  integer, intent(in) :: ievent  ! pre/post sawtooth hint

  real :: ztime_r4,ztime_r4_save

  ztime_r4=ztime
  ztime_r4_save=ztime

  call trx_chk_saw(ztime_r4,ievent)

  if(ztime_r4.ne.ztime_r4_save) then
     ztime = ztime_r4  ! time was adjusted
  endif

end subroutine trx_chk_saw_r8

subroutine trx_chk_saw(ztime,ievent)

  ! check time -- prevent a time that falls between a pre-sawtooth record
  ! and a post-sawtooth record

  ! if ievent=1: favor pre-sawtooth record
  ! if ievent=2: favor post-sawtooth record
  ! otherwise take the nearest record

  use trx_module
  implicit NONE

  real, intent(inout) :: ztime  ! data time (possibly to be adjusted)
  integer, intent(in) :: ievent ! event code (see comments)

  !------------------------------------------
  real, parameter :: zrtol = 2.0e-7 
  real, parameter :: zdtmin= 2.0e-6

  real :: zdtol,zt1,zt2
  integer :: it,jt,ikmax,ikmin
  !------------------------------------------

  if(ztime.lt.time_sc(1)) ztime=time_sc(1)
  if(ztime.gt.time_sc(nsctime)) ztime=time_sc(nsctime)

  do it=1,nsctime-1
     if((time_sc(it).le.ztime).and.(ztime.le.time_sc(it+1))) then
        jt=it
        exit
     endif
  enddo

  ikmin=min(kevent(jt),kevent(jt+1))
  ikmax=max(kevent(jt),kevent(jt+1))
  if(ikmax.eq.0) return

  ! time is in, or adjacent to, a sawtooth interval

  if((ievent.lt.1).or.(ievent.gt.2)) then
     ! no hint
     if(ikmin.eq.0) return  ! not actually in the interval

     ! choose closer time at either end of event transition interval
     if((ztime-time_sc(jt)).le.(time_sc(jt+1)-ztime)) then
        ztime=time_sc(jt)
     else
        ztime=time_sc(jt+1)
     endif
     return
  endif

  ! OK event hint is present; find event range

  zdtol=max(zdtmin,zrtol*max(abs(time_sc(jt)),abs(time_sc(jt+1))))

  if(ikmin.gt.0) then
     ! inside interval
     if(ievent.eq.1) then
        ztime=time_sc(jt)
     else
        ztime=time_sc(jt+1)
     endif
     return
  else if(kevent(jt+1).eq.1) then
     if(ievent.eq.1) return
     !  hint: after sawtooth; accepted only if within zdtol of the interval
     zt1=time_sc(jt+1)-zdtol
     if(ztime.le.zt1) then
        return
     else
        ztime=time_sc(jt+2)
        return
     endif
  else
     if(ievent.eq.2) return
     !  hint: before sawtooth; accept only if within zdtol of the interval
     zt2=time_sc(jt)+zdtol
     if(ztime.ge.zt2) then
        return
     else
        ztime=time_sc(jt-1)
     endif
  endif

end subroutine trx_chk_saw
