subroutine r8_t1scalar(zname,zlabel,zunits,ztime_r8,zdelta_r8,zvalue_r8,ierr)

  ! fetch scalar (REAL*8 interface)
  ! for more info, see t1scalar below...

  use datmgr_mod
  implicit NONE

  character(*), intent(in) :: zname    ! name of scalar function
  character(*), intent(out) :: zlabel  ! function label
  character(*), intent(out) :: zunits  ! function units label

  real*8, intent(in) :: ztime_r8   ! time (seconds) to which to interpolate
  real*8, intent(in) :: zdelta_r8  ! +/- averaging time (seconds)

  real*8, intent(out) :: zvalue_r8 ! the interpolated or averaged value.

  integer, intent(out) :: ierr         ! completion code, 0=OK

  !--------------------------------
  real :: ztime,zdelta,zvalue
  !--------------------------------

  ztime = ztime_r8
  zdelta = zdelta_r8

  call t1scalar(zname,zlabel,zunits,ztime,zdelta,zvalue,ierr)

  zvalue_r8 = zvalue

end subroutine r8_t1scalar

subroutine t1scalar(zname,zlabel,zunits,ztime,zdelta,zvalue,ierr)
!
!  fetch a scalar function of time and interpolate to the indicated
!  time value (ztime), optionally with averaging (+/- zdelta)
!
  use datmgr_mod
  use cplotr_mod
!
  character(*), intent(in) :: zname    ! name of scalar function
  character(*), intent(out) :: zlabel  ! function label
  character(*), intent(out) :: zunits  ! function units label
!
  real, intent(in) :: ztime   ! time (seconds) to which to interpolate
  real, intent(in) :: zdelta  ! +/- averaging time (seconds)
!
  real, intent(out) :: zvalue ! the interpolated or averaged value.
!
  integer, intent(out) :: ierr  ! completion code, 0=OK
!
!-----------------------------
!
!  **local**
!
  character*10 zabbr
!
!-----------------------------
!
  zabbr=zname
  call trcaps(zabbr)
  iadr=ifind_ordr(abt,iordrt,nft,zabbr)
!
  if(iadr.eq.0) then
     call zermsg( &
          ' ?t1scalar:  not a scalar function name:  '//zname)
     ierr=1
     return
  endif
!
  zlabel = labelt(iadr)
  zunits = unitst(iadr)
!
  call dmgfotx(2,ipt,ierr)
  if(ierr.ne.0) then
     ierr=3
     return
  endif
!
  ierr=0
!
  ipf=ipt+(iadr-1)*ntt
!
  ztime1 = ztime-zdelta   ! ztime1==ztime2 possible for small zdelta
  ztime2 = ztime+zdelta
! 
  if((zdelta.le.0.0).or.(ztime+zdelta.le.time3(1)).or. &
       (ztime-zdelta.ge.time3(ntr)) .or. ntr<2 .or. ztime1>=ztime2) then
!
!  straight interpolation
!
     zztime=ztime
     call pltimi_sc(zztime,it1,it2,z1,z2,iexp)
     zvalue = z1*datbuf(ipf+it1-1) + z2*datbuf(ipf+it2-1)
 
  else
!
!  time averaging
!
     zvalue = 0.0
     zwsum = 0.0
     do it=1,ntr-1
        zta=time3(it)
        ztb=time3(it+1)
        ztest1=max(ztime1,zta)
        ztest2=min(ztime2,ztb)
        if(ztest2.gt.ztest1) then
           zztime=0.5*(ztest1+ztest2)
           zdtw=ztest2-ztest1
 
           zwsum=zwsum+zdtw
 
           CALL PLTIMI_SC(ZZTIME,IT1,IT2,Z1,Z2,IEXP)
           zvalue = zvalue + &
                zdtw*(z1*datbuf(ipf+it1-1) + z2*datbuf(ipf+it2-1))
        endif
     enddo

     if (zwsum<=0.) then
        call zermsg( &
             ' ?t1scalar:  error time interpolating on :  '//zname)
        ierr=1
        return
     end if

     zvalue = zvalue/zwsum
 
  endif
!
  return
end subroutine t1scalar
