subroutine tdb_onoff_atime(d,z2char,zthresh,zdtfix,ton,toff,ierr)
  ! for power-channel data (NB,EC,LH,RF) define the on and off times:
  !   ton = min(on times for any channel) (*output*)
  !   toff = max(off times for any channel) (*output*)
  !
  !   z2char => chooses heating channel:
  !     "nb" or "NB" -- neutral beams
  !     "ec" -- ECH/ECCD
  !     "lh" -- Lower Hybrid
  !     "rf" -- ICRF   
  !   all tests of z2char value are case insensitive
  !
  !   zthresh = threshhold:
  !     if positive -- a power, in watts, must be < 20% of the maximum
  !                    power ever occurring on any channel
  !     if negative -- (-1) * a fraction (no units) -- btw -0.0001 and -0.20 --
  !                    power threshold becomes -zthresh * (maximum power
  !                    ever occurring on any channel).
  !
  !   zdtfix -- time to search from first/last powers satisfying thresh-
  !             hold, for an actual 0 or negative power,
  !             or power < threshhold/10
  !
  !   ierr is set only if there is no data or if z2char is unrecognized;
  !        if the threshhold input is out of range it is overridden without
  !        warning
  !
  use trdatbuf_obj
  use tdbsub_uts  ! private
  implicit NONE

  type (trdatbuf) :: d
  character*(*), intent(in) :: z2char  ! channel type NB/LH/EC/RF
  real*8, intent(in) :: zthresh        ! on/off threshhold (see comments)
  real*8, intent(in) :: zdtfix         ! time to search from threshhold
  real*8, intent(out) :: ton,toff      ! on/off times (seconds)
  integer, intent(out) :: ierr         ! completion code; 0=OK

  !-----------------------------------
  real*8, dimension(:,:), allocatable :: zpwrs
  real*8, dimension(:), allocatable :: tonarr,toffarr
  integer :: inchan,intim,ilpwr,iltim,ifound,i
  !-----------------------------------

  ton=ZERO
  toff=ZERO
  ierr=0

  call tdb_pwrset_find(d,z2char,'tdb_onoff_atime',iltim,intim,ilpwr,inchan, &
       ierr)
  if(ierr.gt.0) return
        
  allocate(tonarr(inchan),toffarr(inchan))

  allocate(zpwrs(intim,inchan))
  zpwrs = RESHAPE(d%datbuf(ilpwr:ilpwr+intim*inchan-1), (/ intim, inchan /))
  call tdbsub_onoff(d%datbuf(iltim),intim,zpwrs,inchan, &
       zthresh,zdtfix,tonarr,toffarr)
  deallocate(zpwrs)

  ton=epsinv
  toff=-epsinv
  ifound=0
  do i=1,inchan
     if(tonarr(i).lt.epsinv) then
        ton=min(ton,tonarr(i))
        toff=max(toff,toffarr(i))
        ifound=ifound+1
     endif
  enddo

  if(ifound.eq.0) then
     toff=epsinv*1.001  ! they never come on...
  endif

  deallocate(tonarr,toffarr)
  
end subroutine tdb_onoff_atime

subroutine tdb_onoff_nxtime(d,z2char,zthresh,zdtfix,tprev,tnext,ierr)
  ! for power-channel data (NB,EC,LH,RF) define the "next" on time, after tprev
  !
  ! the "next" on time is a time when power reappears after being zero (or
  ! nearly so) on all channels.  The precise meaning of "nearly so" is in
  ! the threshhold definition.
  !
  !   z2char => chooses heating channel:
  !     "nb" or "NB" -- neutral beams
  !     "ec" -- ECH/ECCD
  !     "lh" -- Lower Hybrid
  !     "rf" -- ICRF   
  !   all tests of z2char value are case insensitive
  !
  !   zthresh = threshhold:
  !     if positive -- a power, in watts, must be < 20% of the maximum
  !                    power ever occurring on any channel
  !     if negative -- (-1) * a fraction (no units) -- btw -0.0001 and -0.20 --
  !                    power threshold becomes -zthresh * (maximum power
  !                    ever occurring on any channel).
  !
  !   zdtfix -- time to search from first/last powers satisfying thresh-
  !             hold, for an actual 0 or negative power, 
  !             or power < threshhold/10
  !
  !   ierr is set only if there is no data or if z2char is unrecognized.
  !        if the threshhold input is out of range it is overridden without
  !        warning
  !
  use trdatbuf_obj
  use tdbsub_uts  ! private
  implicit NONE

  type (trdatbuf) :: d
  character*(*), intent(in) :: z2char  ! channel type NB/LH/EC/RF
  real*8, intent(in) :: zthresh        ! on/off threshhold (see comments)
  real*8, intent(in) :: zdtfix         ! time to search from threshhold
  real*8, intent(in) :: tprev          ! previous time, after which to search

  real*8, intent(out) :: tnext         ! next time of transition from zero to
  !              non-zero power.  If not found, tnext=EPSINV is returned:
  !              a very large time value.

  integer, intent(out) :: ierr         ! completion code; 0=OK

  !-----------------------------------
  integer :: inchan,intim,ilpwr,iltim
  !-----------------------------------

  tnext=EPSINV
  ierr=0

  call tdb_pwrset_find(d,z2char,'tdb_onoff_nxtime',iltim,intim,ilpwr,inchan, &
       ierr)
  if(ierr.gt.0) return
        
  call tdbsub_next_on(d%datbuf(iltim),intim,d%datbuf(ilpwr),inchan, &
       zthresh,zdtfix,tprev,tnext)

end subroutine tdb_onoff_nxtime

subroutine tdb_pwrset_find(d,z2char,zsubr,iltim,intim,ilpwr,inchan,ierr)
  !
  !  find power vs. (time,channel #) dataset
  !
  !   z2char => chooses heating channel:
  !     "nb" or "NB" -- neutral beams
  !     "ec" -- ECH/ECCD
  !     "lh" -- Lower Hybrid
  !     "rf" -- ICRF   
  !   all tests of z2char value are case insensitive
  !

  use trdatbuf_obj
  use tdbsub_uts  ! private
!  use trcom, only: imas_ref_shot
  implicit NONE

  type (trdatbuf) :: d
  character*(*), intent(in) :: z2char  ! channel type NB/LH/EC/RF
  character*(*), intent(in) :: zsubr   ! caller ID string (for error message)
  integer, intent(out) :: iltim        ! address of time data
  integer, intent(out) :: intim        ! number of time points
  integer, intent(out) :: ilpwr        ! address of power data
  integer, intent(out) :: inchan       ! number of power channels
  integer, intent(out) :: ierr         ! status code returned: 0=normal
  !
  !  note ierr=1 is returned if inchan=0, but, this could be a normal situation
  !  indicating e.g. a simulation where there is no heating of the type 
  !  requested
  !
  !-----------------------------------------
  integer :: lunmsg_tdb
  !-----------------------------------------
  !  local variables:
  character*3 :: ztest3
  !-----------------------------------------

  iltim=0
  intim=0
  ilpwr=0
  inchan=0
  ierr=0

  ztest3=z2char(1:min(3,len(z2char)))
  call uupper(ztest3)

  !  look for antenna or NB set...

  !  expect channel based data (NB powers, RF antennas etc)

  if(ztest3.eq.'NB') then
     inchan=d%nbdata
     if(inchan.eq.0) then
        ierr=1
     else
        ilpwr=d%lpwrnb
        iltim=d%ltimnb
        intim=d%ntimnb
     endif
  else if(ztest3.eq.'EC') then
     inchan=d%nantech_d
     if(inchan.eq.0) then
        ierr=1
     else
        ilpwr=d%lpwrec
        ! having issues introducing new flag in datbuf, so we check the size
        ! until that is fixed
        if (d%lfecq.le.0) then
!        if (imas_ref_shot.le.0) then
           iltim=d%ltimec
           intim=d%ntimec
        else 
           ! data come from IDS, interpolated over time1
           iltim=d%ltime1
           intim=d%ntime1
        endif
     endif
  else if(ztest3.eq.'LH') then
     inchan=d%nantlh_d
     if(inchan.eq.0) then
        ierr=1
     else
        ilpwr=d%lpwrlh
        iltim=d%ltimlh
        intim=d%ntimlh
     endif
  else if(ztest3.eq.'RF') then
     inchan=d%nantich_d
     if(inchan.eq.0) then
        ierr=1
     else
        ilpwr=d%lpwrrf
        iltim=d%ltimrf
        intim=d%ntimrf
     endif
  else
     ierr=2
     write(lunmsg_tdb(0),*) &
          ' ? trdatbuf_lib/'//trim(zsubr)//': unrecognized ', &
          'channel identifier: "',trim(z2char),'".'
  endif

end subroutine tdb_pwrset_find

subroutine tdb_onoff_atime1(d,ntri,ztris,zthresh,zdtfix,ton,toff,ierr)
  ! for a set of time series data Pj(t) define the on and off times:
  !   ton = min(on times for any channel) (*output*)
  !   toff = max(off times for any channel) (*output*)
  !
  !   ztris(1:ntri) -- trigraphs of scalar functions
  !
  !   zthresh = threshhold:
  !     if positive -- a power, in watts, must be < 20% of the maximum
  !                    power ever occurring on any channel
  !     if negative -- (-1) * a fraction (no units) -- btw -0.0001 and -0.20 --
  !                    power threshold becomes -zthresh * (maximum power
  !                    ever occurring on any channel).
  !
  !   zdtfix -- time to search from first/last powers satisfying thresh-
  !             hold, for an actual 0 or negative power...
  !
  !   ierr is set only if there is no data or if z2char is unrecognized;
  !   if the zthresh limit has to be adjusted to conform to rules, a 
  !   warning message is written but ierr is not set.
  !
  use trdatbuf_obj
  use tdbsub_uts  ! private
  implicit NONE

  type (trdatbuf) :: d
  integer, intent(in) :: ntri          ! no. of trigraphs
  character*(*), intent(in) :: ztris(ntri)  ! trigraph names
  real*8, intent(in) :: zthresh        ! on/off threshhold (see comments)
  real*8, intent(in) :: zdtfix         ! time to search from threshhold
  real*8, intent(out) :: ton,toff      ! on/off times (seconds)
  integer, intent(out) :: ierr         ! completion code; 0=OK

  !-----------------------------------
  real*8, dimension(:), allocatable :: tonarr,toffarr
  real*8, dimension(:,:), allocatable :: zpwrs
  character*3 :: ztest3
  integer :: inchan,intim,iltim,lunmsg_tdb,ifound,i
  integer :: iadscal
  !-----------------------------------

  ton=epsinv
  toff=1.01*epsinv

  ierr=0

  iltim = d%ltime1
  intim = d%ntime1
  inchan = ntri

  allocate(zpwrs(intim,ntri)); zpwrs=ZERO
  ifound = 0

  do i=1,ntri
     ztest3=ztris(i)
     call uupper(ztest3)

     if(tdb_defined1(d,ztest3)) then

     !  look for scalar function...

        if(tdb_present1(d,ztest3,iadscal)) then

           ifound = ifound + 1
           zpwrs(1:intim,i) = d%datbuf(iadscal:iadscal+intim-1)

        endif

     else

        write(lunmsg_tdb(0),*) ' ?tdb_onoff_atime1: not a scalar function: ', &
             ztest3
        ierr=ierr+1
        
     endif
  enddo

  if((ierr.gt.0).or.(ifound.eq.0)) then
     deallocate(zpwrs)
     return
  endif
        
  allocate(tonarr(inchan),toffarr(inchan))

  call tdbsub_onoff(d%datbuf(iltim),intim,zpwrs,inchan, &
       zthresh,zdtfix,tonarr,toffarr)
  deallocate(zpwrs)

  ton=epsinv
  toff=-epsinv
  ifound=0
  do i=1,inchan
     if(tonarr(i).lt.epsinv) then
        ton=min(ton,tonarr(i))
        toff=max(toff,toffarr(i))
        ifound=ifound+1
     endif
  enddo

  if(ifound.eq.0) then
     toff=epsinv*1.01  ! they never come on...
  endif

  deallocate(tonarr,toffarr)
  
end subroutine tdb_onoff_atime1

subroutine tdb_onoff_atime2(d,ntri,ztris,zthresh,zdtfix,ton,toff,ierr)
  ! for a set of profiles Pj(t,x) define the on and off times:
  !   ton = min(on times for any profile) (*output*)
  !   toff = max(off times for any profile) (*output*)
  !
  !   ztris(1:ntri) -- trigraphs of scalar functions
  !
  !   zthresh = threshhold:
  !     if non-negative -- same as zthresh=-0.05 here
  !        (total power not usable, we don't know it here)
  !
  !     if negative -- (-1) * a fraction (no units) -- btw -0.0001 and -0.20 --
  !                    power threshold becomes -zthresh * (maximum power
  !                    ever occurring on any channel).
  !
  !   zdtfix -- time to search from first/last powers satisfying thresh-
  !             hold, for an actual 0 or negative power...
  !
  !   ierr is set only if the trigraphs are unrecognized.
  !
  use trdatbuf_obj
  use tdbsub_uts  ! private
  implicit NONE

  type (trdatbuf) :: d
  integer, intent(in) :: ntri          ! no. of trigraphs
  character*(*), intent(in) :: ztris(ntri)  ! trigraph names
  real*8, intent(in) :: zthresh        ! on/off threshhold (see comments)
  real*8, intent(in) :: zdtfix         ! time to search from threshhold
  real*8, intent(out) :: ton,toff      ! on/off times (seconds)
  integer, intent(out) :: ierr         ! completion code; 0=OK

  !-----------------------------------
  real*8, dimension(:), allocatable :: tonarr,toffarr
  real*8, dimension(:,:),  allocatable :: zpwrs
  real*8 :: zthreshi
  character*3 :: ztest3
  integer :: inchan,intim,iltim,lunmsg_tdb,ifound,i,it,ix,inx
  integer :: iadprof,iadscal,iadxpro,inxpro,iadxt
  !-----------------------------------

  ton=epsinv
  toff=1.01*epsinv
  ierr=0

  iltim=d%ltime2
  intim=d%ntime2
  inchan=ntri

  allocate(zpwrs(intim,inchan)); zpwrs = ZERO
  ifound = 0

  do i=1,ntri
     ztest3=ztris(i)
     call uupper(ztest3)

     if(tdb_defined2(d,ztest3)) then

        !  look for profile function...

        if(tdb_present2(d,ztest3,iadprof,iadxpro,inxpro)) then

           ifound = ifound+1

           if(inxpro.eq.0) then
              ierr=3
              write(lunmsg_tdb(0),*) &
                   ' ?? tdb_onoff_atime: cannot infer on/off time from: '//ztest3
              write(lunmsg_tdb(0),*) '    profile type not supported.'
           else

              do ix=1,inxpro
                 do it=1,intim

                    iadxt=iadprof+(ix-1)*intim + it - 1
                    zpwrs(it,i) = zpwrs(it,i) + d%datbuf(iadxt)

                 enddo
              enddo
           endif

        endif

     else

        write(lunmsg_tdb(0),*) ' ?tdb_onoff_atime2: not a profile function: ',&
             ztest3
        ierr=1

     endif
  enddo

  if((ierr.gt.0).or.(ifound.eq.0)) then
     deallocate(zpwrs)
     return
  endif
        
  allocate(tonarr(inchan),toffarr(inchan))

  zthreshi=zthresh
  if(zthreshi.ge.ZERO) zthreshi=-0.05d0

  call tdbsub_onoff(d%datbuf(iltim),intim,zpwrs,inchan, &
       zthreshi,zdtfix,tonarr,toffarr)
  deallocate(zpwrs)

  ton=epsinv
  toff=-epsinv
  ifound=0
  do i=1,inchan
     if(tonarr(i).lt.epsinv) then
        ton=min(ton,tonarr(i))
        toff=max(toff,toffarr(i))
        ifound=ifound+1
     endif
  enddo

  if(ifound.eq.0) then
     toff=epsinv*1.01  ! they never come on...
  endif

  deallocate(tonarr,toffarr)
  
end subroutine tdb_onoff_atime2

subroutine tdb_onoff_times(d,z2char,zthresh,zdtfix,tonarr,toffarr,ierr)
  ! for power-channel data (NB,EC,LH,RF) define the on and off times:
  !   tonarr(i) = on time for channel #i (i.e. beam or antenna) (*output*)
  !   toffarr(i) = off time for channel #i (*output*)
  !
  !   z2char => chooses heating channel:
  !     "nb" or "NB" -- neutral beams
  !     "ec" -- ECH/ECCD
  !     "lh" -- Lower Hybrid
  !     "rf" -- ICRF   
  !   all tests of z2char value are case insensitive
  ! 
  !   zthresh = threshhold:
  !     if positive -- a power, in watts, must be < 20% of the maximum
  !                    power ever occurring on any channel
  !     if negative -- (-1) * a fraction (no units) -- btw -0.0001 and -0.20 --
  !                    power threshold becomes -zthresh * (maximum power
  !                    ever occurring on any channel).
  !
  !   zdtfix -- time to search from first/last powers satisfying thresh-
  !             hold, for an actual 0 or negative power...
  !
  !   ierr is set only if there is no data or if z2char is unrecognized;
  !   if the zthresh limit has to be adjusted to conform to rules, a 
  !   warning message is written but ierr is not set.
  !
  use trdatbuf_obj
  use tdbsub_uts  ! private
  implicit NONE

  type (trdatbuf) :: d
  character*(*), intent(in) :: z2char  ! channel type NB/LH/EC/RF
  real*8, intent(in) :: zthresh        ! on/off threshhold (see comments)
  real*8, intent(in) :: zdtfix         ! time to search from threshhold
  real*8, intent(out), dimension(:) :: tonarr,toffarr  ! on/off times (seconds)
  integer, intent(out) :: ierr         ! completion code; 0=OK

  !-----------------------------------
  real*8, dimension(:,:), allocatable :: zpwrs
  integer :: inchan,intim,ilpwr,iltim,lunmsg_tdb,ifound,i
  !-----------------------------------

  ierr=0

  call tdb_pwrset_find(d,z2char,'tdb_onoff_times',iltim,intim,ilpwr,inchan, &
       ierr)

  if(size(tonarr).lt.inchan) then
     write(lunmsg_tdb(0),*) ' ? trdatbuf_lib/tdb_onoff_times: "tonarr" array',&
          ' size too small:'
     write(lunmsg_tdb(0),*) '   need: ',inchan
     write(lunmsg_tdb(0),*) '    got: ',size(tonarr)
     ierr=3
  endif

  if(size(toffarr).lt.inchan) then
     write(lunmsg_tdb(0),*) ' ? trdatbuf_lib/tdb_onoff_times: "toffarr" array',&
          ' size too small:'
     write(lunmsg_tdb(0),*) '   need: ',inchan
     write(lunmsg_tdb(0),*) '    got: ',size(toffarr)
     ierr=4
  endif

  if(ierr.gt.0) return

  allocate(zpwrs(intim,inchan))
  zpwrs = RESHAPE(d%datbuf(ilpwr:ilpwr+intim*inchan-1), (/ intim, inchan /))
  call tdbsub_onoff(d%datbuf(iltim),intim,zpwrs,inchan, &
       zthresh,zdtfix,tonarr(1:inchan),toffarr(1:inchan))
  deallocate(zpwrs)
end subroutine tdb_onoff_times

subroutine tdb_onoff_zscal(d,ztri,zton,ztoff,ierr)

  !  zero out a scalar function outside a given time range
  use trdatbuf_obj
  use tdbsub_uts  ! private
  implicit NONE

  type (trdatbuf) :: d
  character*(*), intent(in) :: ztri  ! trigraph of data to modify
  real*8, intent(in) :: zton,ztoff   ! on/off times
  integer, intent(out) :: ierr       ! 0 if OK, set if ztri invalid...

  !------------------------
  integer :: iaddr,iltim,intim,it
  !------------------------

  if(tdb_present1(d,ztri,iaddr)) then
     iltim=d%ltime1
     intim=d%ntime1
     do it=1,intim
        if((d%datbuf(iltim+it-1).lt.zton).or. &
             (d%datbuf(iltim+it-1).gt.ztoff)) then
           d%datbuf(iaddr+it-1) = ZERO
        endif
     enddo

     ierr = 0

  else

     ierr = 1

  endif
end subroutine tdb_onoff_zscal

subroutine tdb_onoff_zprof(d,ztri,zton,ztoff,ierr)

  !  zero out a profile function outside a given time range
  use trdatbuf_obj
  use tdbsub_uts  ! private
  implicit NONE

  type (trdatbuf) :: d
  character*(*), intent(in) :: ztri  ! trigraph of data to modify
  real*8, intent(in) :: zton,ztoff   ! on/off times
  integer, intent(out) :: ierr       ! 0 if OK, set if ztri invalid...

  !------------------------
  integer :: iaddr,iaddx,inumx,iltim,intim,it,ix
  !------------------------

  if(tdb_present2(d,ztri,iaddr,iaddx,inumx)) then
     if(inumx.eq.0) then
        ierr = 2
     else
        iltim=d%ltime2
        intim=d%ntime2
        do it=1,intim
           if((d%datbuf(iltim+it-1).lt.zton).or. &
                (d%datbuf(iltim+it-1).gt.ztoff)) then
              do ix=1,inumx
                 d%datbuf(iaddr+(ix-1)*intim+it-1) = ZERO
              enddo
           endif
        enddo
        ierr = 0
     endif

  else

     ierr = 1

  endif
end subroutine tdb_onoff_zprof

