subroutine tdbsub_onoff(time,intime,pdata,inchan,zthresh,zdtfix,ton,toff)

  use tdbsub_uts
  implicit NONE

  !  analyze power channel data for on/off times
  !  see tdb_onoff_times.f90 comments for definition of "zthresh" & "zdtfix"

  integer, intent(in) :: intime  ! # of times
  real*8, intent(in) :: time(intime)  ! the time values
  integer, intent(in) :: inchan  ! # of channels (i.e beams or RF antennas)
  real*8, intent(in) :: pdata(intime,inchan)

  real*8, intent(in) :: zthresh  ! power threshhold, watts or -fraction
  real*8, intent(in) :: zdtfix   ! time to search around threshhold points

  real*8, intent(out) :: ton(inchan)  ! inferred on times
  real*8, intent(out) :: toff(inchan) ! inferred off times

  !----------------------------------------
  real*8 :: zthresha
  integer :: ichan,it,it1,it2,ifound
  real*8 :: ztfound,zpfound
  !----------------------------------------

  call tdbsub_onoff_thresh(intime,pdata,inchan,zthresh,zthresha)

  do ichan=1,inchan
     !  find first time with power above threshhold
     ifound=0
     do it=1,intime
        if(pdata(it,ichan).gt.zthresha) then
           ifound=it
           exit
        endif
     enddo

     if(ifound.eq.0) then
        ! if nothing is found, set on/off times to very large numbers
        ton(ichan)=EPSINV
        toff(ichan)=EPSINV
        cycle
     endif

     ztfound=time(ifound)
     zpfound=zthresha
     do it=ifound,1,-1
        it1=min(ifound,(it+1))
        if((ztfound-time(it)).ge.zdtfix) exit
        if(pdata(it,ichan).lt.zpfound*0.1d0) then
           it1=it+1
           exit
        endif
     enddo

     if(it1.eq.ifound) then
        ton(ichan)=time(max(1,(it1-1)))
     else
        ton(ichan)=HALF*(time(it1)+time(max(1,it1-1)))
     endif

     !  find last time with power above threshhold -- something exists
     !  or the loop would have exited above...

     ifound=0
     do it=intime,1,-1
        if(pdata(it,ichan).gt.zthresha) then
           ifound=it
           exit
        endif
     enddo

     ztfound=time(ifound)
     zpfound=pdata(ifound,ichan)
     do it=ifound,intime
        it2=max(ifound,(it-1))
        if((time(it)-ztfound).ge.zdtfix) exit
        if(pdata(it,ichan).lt.zpfound*0.1d0) then
           it2=it-1
           exit
        endif
     enddo

     if(it2.eq.ifound) then
        toff(ichan)=time(min(intime,it2+1))
     else
        toff(ichan)=HALF*(time(it2)+time(min(intime,it2+1)))
     endif
  enddo

end subroutine tdbsub_onoff

subroutine tdbsub_next_on(time,intime,pdata,inchan,zthresh,zdtfix,tprev,tnext)

  use tdbsub_uts
  implicit NONE

  !  find the first "on time" for power data, following time tprev
  !  This is defined as a time when one or more channels are "on", immediately
  !  following a time (.gt.tprev) when all are off.

  !  mod DMC Feb 2010 -- bugfix -- for all to count as off, they must remain
  !  off for at least 2 consecutive time points-- otherwise false on/off
  !  events were detectable with very close times, causing time step
  !  controller problems...

  !  see tdb_onoff_times.f90 comments for definition of "zthresh" & "zdtfix"

  !  NOTE: it is assumed, without checking, that the times are in strict
  !  ascending order.

  integer, intent(in) :: intime  ! # of times
  real*8, intent(in) :: time(intime)  ! the time values
  integer, intent(in) :: inchan  ! # of channels (i.e beams or RF antennas)
  real*8, intent(in) :: pdata(intime,inchan)

  real*8, intent(in) :: zthresh  ! power threshhold, watts or -fraction
  real*8, intent(in) :: zdtfix   ! time to search around threshhold points

  real*8, intent(in) :: tprev    ! prior time; return tnext > tprev
  real*8, intent(out) :: tnext   ! time of next transition from zero to finite
  !  power...

  !----------------------------------------
  real*8 :: zthresha
  integer :: ichan,it,it0,it1,it2,ifound,icount
  real*8 :: ztfound,zpfound,zptest
  integer :: ichan_with_power(inchan)
  !----------------------------------------

  call tdbsub_onoff_thresh(intime,pdata,inchan,zthresh,zthresha)

  if(tprev.ge.time(intime)) then
     ! no more times
     tnext = EPSINV
     return
  endif

  ! we now know that a time t > tprev exists; find the first one...

  if(tprev.lt.time(1)) then
     it0=1
     it1=1  ! just set up to find the first ontime...

  else
     ! (tprev.eq.time(1)) lands here...

     do it=1,intime
        if(time(it).gt.tprev) then
           it0=it
           exit
        endif
     enddo

     ! find the first time, at/after it0, where all the powers are below
     ! the threshhold for 2 consecutive time points

     it1 = -1
     do it=it0,intime-1
        ifound=0
        do ichan=1,inchan
           if(max(pdata(it,ichan),pdata(it+1,ichan)).gt.zthresha) then
              ifound=it
              exit
           endif
        enddo

        if(ifound.eq.0) then
           it1=it
           exit
        endif
     enddo

     if(it1.eq.-1) then
        ! no times with near zero power, after tprev
        tnext = EPSINV
        return
     endif

  endif

  ! find the first time, after it1, where one or more powers are above
  ! the threshhold

  it2 = -1
  ichan_with_power = 0

  do it=it1,intime
     ifound=0
     do ichan=1,inchan
        if(pdata(it,ichan).gt.zthresha) then
           ifound=it
           ichan_with_power(ichan)=1
        endif
     enddo

     if(ifound.gt.0) then
        it2=it
        exit
     endif
  enddo

  if(it2.eq.-1) then
     ! no times with power above threshhold, after zero after tprev
     tnext = EPSINV
     return
  endif

  ifound=it2
  ztfound=time(ifound)
  zpfound=zthresha

  do it=ifound,1,-1
     it1=min(ifound,(it+1))
     if((ztfound-time(it)).ge.zdtfix) exit

     icount=0
     do ichan=1,inchan
        if(ichan_with_power(ichan).eq.1) then
           if(pdata(it,ichan).ge.zpfound*0.1d0) then
              icount=icount+1
           else
              ichan_with_power(ichan)=0
           endif
        endif
     enddo

     if(icount.eq.0) then
        it1=it+1
        exit
     endif
  enddo

  if(it1.eq.ifound) then
     tnext=time(max(1,(it1-1)))
  else
     tnext=HALF*(time(it1)+time(max(1,it1-1)))
  endif

end subroutine tdbsub_next_on

subroutine tdbsub_onoff_thresh(intime,pdata,inchan,zthresh,zthresha)

  use tdbsub_uts
  implicit NONE

  ! get the absolute power threshold
  ! (see tdb_onoff_times.f90 comments...)

  integer, intent(in) :: intime   ! no. of time points in channel data
  integer, intent(in) :: inchan   ! no. of channels

  real*8, intent(in) :: pdata(intime,inchan)  ! channel data (power)

  real*8, intent(in) :: zthresh   ! threshhold control

  real*8, intent(out) :: zthresha ! absolute threshold (in power units)

  !----------------------------------------
  real*8 :: zmax
  real*8, parameter :: zmaxf=0.20d0
  real*8, parameter :: zminf=0.0001d0
  integer :: ichan,it,it1,it2,ifound
  real*8 :: ztfound,zpfound
  !----------------------------------------

  zmax=maxval(pdata)
  zmax=max(ONE,zmax)  ! at least one watt, please...

  if(zthresh.ge.ZERO) then
     ! threshhold in watts, but constrain according to data extrema
     zthresha = max(zminf*zmax,min(zmaxf*zmax,zthresh))
  else
     ! fraction
     zthresha=max(zminf,min(zmaxf,-zthresh))
     zthresha=zthresha*zmax
  endif

end subroutine tdbsub_onoff_thresh

