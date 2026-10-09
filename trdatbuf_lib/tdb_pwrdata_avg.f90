subroutine tdb_pwrdata_avg(d,zpwr,ierr)
  use trdatbuf_obj
  use trdatbuf_aux
  implicit NONE

  !  get beam parameters averaged over a time step
  !  Mods
  !     dmc Aug 2007: if a time step of zero length is provided, 
  !                   impose instead a very small but finite time step).
  !
  !     08Aug2011     (to rev 15166 )    jim.conboy@ccfe.ac.uk
  !                   Revise error msgs, identify item requested
  !

  type (trdatbuf) :: d   ! trdat data object
  type (pwrget) :: zpwr  ! specification of desired parameters
  integer, intent(out) :: ierr   ! completion code (0=OK)

  integer :: lunmsg_tdb
  integer :: ipwr,iparam,inb,int,iadr,iadrp,ilt,ib,inbx
  real*8 :: zt1,zt2,zfdt,zpfdt,zpdt,ztdat1,ztdat2
  real*8, parameter :: ZERO = 0.0d0
  real*8, parameter :: ZEPS = 1.0d-12
  real*8 :: tdbsub_i1,tdbsub_i2

  real*8 :: dtmin_fac,dtmin,zdelt

  character*3 zitem

  logical :: shortn,datlim
  character(len=*), parameter   :: cr="trdatbuf_lib/tdb_pwrdata_avg"
  !------------------------------------------------------
  !  error checks...
  ierr=0    ! think positive

  zpwr%zparam = ZERO

  !  zpwr is expected to have been initialized...
  if(zpwr%nguard.ne.123456789) then
     ierr=1
     write(lunmsg_tdb(0),*) &
          ' ?'//cr//': (pwrget) object was not initialized.'
  endif

  zitem=zpwr%item
  call uupper(zitem)

  shortn=.FALSE.
  datlim=.FALSE.

  if(zitem.eq.'PWR') then
     ipwr = d%lpwrnb
     inbx = d%nbdata
     iparam = d%lpwrnb
     ilt = d%ltimnb
     int = d%ntimnb

  else if(zitem.eq.'VLT') then
     ipwr = d%lpwrnb
     inbx = d%nbdata
     iparam = d%lvltnb
     ilt = d%ltimnb
     int = d%ntimnb

  else if(zitem.eq.'FUL') then
     ipwr = d%lpwrnb
     inbx = d%nbdata
     iparam = d%lfulnb
     ilt = d%ltimnb
     int = d%ntimnb

  else if(zitem.eq.'HLF') then
     ipwr = d%lpwrnb
     inbx = d%nbdata
     iparam = d%lhlfnb
     ilt = d%ltimnb
     int = d%ntimnb

  else if(zitem.eq.'RFP') then
     ipwr = d%lpwrrf
     inbx = d%nantich_d
     iparam = d%lpwrrf
     ilt = d%ltimrf
     int = d%ntimrf

  else if(zitem.eq.'RFF') then
     ipwr = d%lfrqrff    ! timebase difference; cannot use lpwrrf
     shortn = .TRUE.     ! constrain time average btw on & off times
     inbx = d%nantich_d
     iparam = d%lfrqrff
     ilt = d%ltimrff
     int = d%ntimrff

  else if(zitem.eq.'ECP') then
     ipwr = d%lpwrec
     inbx = d%nantech_d
     iparam = d%lpwrec
     ilt = d%ltimec
     int = d%ntimec

  else if(zitem.eq.'ECA') then
     ipwr = d%lfeca      ! timebase difference; cannot use lpwrec
     datlim = .TRUE.     ! constrain time average btw data start/stop times
     inbx = d%nantech_d
     iparam = d%lfeca
     ilt = d%ltimeca
     int = d%ntimeca

  else if(zitem.eq.'ECB') then
     ipwr = d%lfecb      ! timebase difference; cannot use lpwrec
     datlim = .TRUE.     ! constrain time average btw data start/stop times
     inbx = d%nantech_d
     iparam = d%lfecb
     ilt = d%ltimecb
     int = d%ntimecb

  else if(zitem.eq.'LHP') then
     ipwr = d%lpwrlh
     inbx = d%nantlh_d
     iparam = d%lpwrlh
     ilt = d%ltimlh
     int = d%ntimlh

  else
     ierr=3
     write(lunmsg_tdb(0),*) &
          ' ?'//cr//': unrecognized data item: "', zitem,'"'
  endif

  if(ipwr.eq.0) then
     ierr=2
     write(lunmsg_tdb(0),*) &
          ' ?'//cr//': called with no power data available -',zitem
  endif

  if(zpwr%ztime1.gt.zpwr%ztime2) then
     write(lunmsg_tdb(0),*) &
          ' ?'//cr//': averaging times out of order -',zitem
     write(lunmsg_tdb(0),*) &
          '  zpwr%ztime1 = ',zpwr%ztime1,'   zpwr%ztime2 = ',zpwr%ztime2
     if( zpwr%ztime1==0. .and.  zpwr%ztime2== -1. )  &
         write(lunmsg_tdb(0),*) &
           zitem,' parameters were  never set (?) '
     ierr=4
  endif

  if(zpwr%nbeam .ne. inbx) then
     write(lunmsg_tdb(0),*) ' ?',cr,': ',zitem, &
          '-  inconsistent #beams or #antennas in zpwr & d objects:'
     write(lunmsg_tdb(0),*) &
          '  in d = ',inbx,'   in zpwr = ',zpwr%nbeam
     ierr=5
  endif

  inb = zpwr%nbeam

  if(zpwr%tbon(1).gt.zpwr%tboff(1)) then
     write(lunmsg_tdb(0),*) &
          ' ?'//cr//': ',zitem,' - on/off times out of order:'
     write(lunmsg_tdb(0),*) &
          '  zpwr%tbon(1) = ',zpwr%tbon(1),'   zpwr%tboff(1) = ',zpwr%tboff(1)
     if( zpwr%tbon(1) == 0. .and. zpwr%tboff(1)== -1. )  &
        write(lunmsg_tdb(0),*) &
           zitem,' parameters were  never set (?) '
     ierr=6
  endif

  if(ierr.gt.0) then
     return
  endif

  !-----------
  !  OK...

  do ib=1,inb
     if(datlim) then
        ztdat1 = d%datbuf(ilt)
        ztdat2 = d%datbuf(ilt+int-1)
        zt1=max(ztdat1,min(ztdat2,zpwr%ztime1))
        zt2=max(ztdat1,min(ztdat2,zpwr%ztime2))
     else
        zt1=max(zpwr%ztime1,zpwr%tbon(ib))
        zt2=min(zpwr%ztime2,zpwr%tboff(ib))
     endif

     if(zt2.lt.zt1) then
        zpwr%zparam(ib)=ZERO

     else if( ((zt2.lt.zpwr%tbon(ib)).or.(zt1.gt.zpwr%tboff(ib))) .AND. &
          (.not.datlim) ) then
        zpwr%zparam(ib)=ZERO

     else
        dtmin_fac=zeps
        dtmin=max(zeps,zeps*max(abs(zt1),abs(zt2)))
        if(zt2.lt.zt1+dtmin) then
           zt2=zt1+dtmin
        endif

        if(zpwr%pweight) then
           iadrp = ipwr + (ib-1)*int
           iadr = iparam + (ib-1)*int
           zpdt = tdbsub_i1(zt1,zt2,d%datbuf(ilt),int,d%datbuf(iadrp))
           zpfdt = tdbsub_i2(zt1,zt2,d%datbuf(ilt),int, &
                d%datbuf(iadrp),d%datbuf(iadr))
           if(zpdt.le.ZERO) then
              zpwr%zparam(ib)=ZERO
           else
              zpwr%zparam(ib)=zpfdt/zpdt
           endif
        else
           iadr = iparam + (ib-1)*int
           zfdt = tdbsub_i1(zt1,zt2,d%datbuf(ilt),int,d%datbuf(iadr))
           if(shortn.or.datlim) then
              zdelt=(zt2-zt1)
           else
              zdelt=max(dtmin,(zpwr%ztime2-zpwr%ztime1))
           endif
           zpwr%zparam(ib)=zfdt/zdelt
        endif
     endif
  enddo

end subroutine tdb_pwrdata_avg
