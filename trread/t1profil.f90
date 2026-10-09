subroutine r8_t1profil(zname,zlabel,zunits,ztime_r8,zdelta_r8, &
     istype,zprofil_r8,nmax,ngot,ierr)

  ! fetch profile -- REAL*8 interface; see t1profil below.

  use datmgr_mod
  implicit NONE

  character(*), intent(in) :: zname    ! name of profile function
  character(*), intent(out) :: zlabel  ! function label
  character(*), intent(out) :: zunits  ! function units label

  real*8, intent(in) :: ztime_r8   ! time (seconds) to which to interpolate
  real*8, intent(in) :: zdelta_r8  ! +/- averaging time (seconds)

  integer, intent(out) :: istype   ! type of x-axis associated with profile
  integer, intent(in) :: nmax      ! dimension of zprofil array
  real*8, dimension(nmax) :: zprofil_r8 ! the output profile (1:ngot written)
  integer, intent(out) :: ngot     ! actual size of profile returned

  integer, intent(out) :: ierr     ! completion code, 0=OK

  !-------------------------------------------------------
  real :: ztime,zdelta
  real, dimension(:), allocatable :: zprofil
  integer :: ii
  !-------------------------------------------------------

  ztime = ztime_r8
  zdelta = zdelta_r8

  allocate(zprofil(nmax))

  call t1profil(zname,zlabel,zunits,ztime,zdelta, &
     istype,zprofil,nmax,ngot,ierr)

  do ii=1,ngot
     zprofil_r8(ii) = zprofil(ii)
  enddo

end subroutine r8_t1profil


subroutine t1profil(zname,zlabel,zunits,ztime,zdelta, &
     istype,zprofil,nmax,ngot,ierr)

  use datmgr_mod
  use cplotr_mod

!  fetch a profile function of time and interpolate to the indicated
!  time value (ztime), optionally with averaging (+/- zdelta)
!
  character(*), intent(in) :: zname    ! name of profile function
  character(*), intent(out) :: zlabel  ! function label
  character(*), intent(out) :: zunits  ! function units label
!
  real, intent(in) :: ztime   ! time (seconds) to which to interpolate
  real, intent(in) :: zdelta  ! +/- averaging time (seconds)
!
  integer, intent(out) :: istype  ! type of x-axis associated with profile
  integer, intent(in) :: nmax   ! dimension of zprofil array
  real, dimension(nmax) :: zprofil   ! the output profile (1:ngot written)
  integer, intent(out) :: ngot  ! actual size of profile returned
!
  integer, intent(out) :: ierr  ! completion code, 0=OK
!
!-----------------------------
!
!  **local**
!
  character(10) :: zabbr
!
!-----------------------------
!
  zabbr=zname
  call trcaps(zabbr)
  iadr=ifind_ordr(abr,iordrr,nfxt,zabbr)
!
  if(iadr.eq.0) then
     call zermsg( &
          ' ?t1profil:  not a profile function name:  '//zname)
     ierr=1

     istype=-99
     ngot=0

     return
  endif
!
  zlabel = labelr(iadr)
  zunits = unitsr(iadr)
  istype = itypr(iadr)
!
  inx = nzonex(istype)
  if(inx.gt.nmax) then
     call zermsg( &
          ' ?t1profil:  passed "zprofil" array too small.')
     ierr=2

     istype=-99
     ngot=0

     return
  endif
!
  ngot=inx
  ierr=0
!
!  read/find the data
!
  CALL DMGFXT(IADR,IND)
  IPF=LOCD(IND)
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
     call pltimi(zztime,it1,it2,z1,z2,iexp)
     do ix=1,inx
        zprofil(ix)= &
             z1*datbuf(ipf+(it1-1)*inx+ix-1) +z2*datbuf(ipf+(it2-1)*inx+ix-1)
     enddo
 
  else
!
!  time averaging
!
     zprofil(1:ngot)=0.0
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
 
           CALL PLTIMI(ZZTIME,IT1,IT2,Z1,Z2,IEXP)
           do ix=1,inx
              zprofil(ix)= zprofil(ix) + zdtw* &
                   (z1*datbuf(ipf+(it1-1)*inx+ix-1) + &
                    z2*datbuf(ipf+(it2-1)*inx+ix-1))
           enddo
 
        endif
     enddo
 
     if (zwsum<=0.) then
        call zermsg( &
             ' ?t1profil:  error time interpolating on :  '//zname)
        ierr=1
        return
     end if

     zprofil(1:ngot)=zprofil(1:ngot)/zwsum
 
  endif

  return
end subroutine t1profil
