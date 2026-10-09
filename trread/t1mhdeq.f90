subroutine t1mhdeq(ztime,zdelta,nsmax,ntheta,nsgot,units_option, &
     theta, &
     Rarr,Zarr, &
     rho,psi,pmhd,qmhd,gmhd, &
     tflux, pcur, ierr )
!
  implicit none
!
!  read in the TRANSP MHD equilibrium at a particular time
!  or averaged over a range of times
!
! in seconds:
!
  real, intent(in) :: ztime           ! time at which equilibrium is wanted
  real, intent(in) :: zdelta          ! +/- zdelta for time average, or 0.0
!
! dimensioning info:
!
  integer, intent(in) :: nsmax        ! max no. of surfaces, counting axis
  integer, intent(in) :: ntheta       ! no. of theta points in Rarr,Zarr
  integer, intent(out) :: nsgot       ! actual no. of surfaces counting axis
!
! output physical units option
!
  character(*), intent(in) :: units_option
!                             "TRANSP" for TRANSP traditional units
!                             "MKS"    for MKS units:
!              TRANSP  MKS
!  Rarr,Zarr     cm     m
!   rho                      (dimensionless)
!   psi         Webers/rad   (same for both)
!   Pmhd        Pascals      (same for both)
!   Qmhd                     (dimensionless)
!   Gmhd        T*cm   T*m
!   Tflux       Webers       (same for both)
!   Pcura       Amps         (same for both)
!
  real, intent(in), dimension(ntheta) :: theta         ! theta grid, input
!---------------------------
!  equilibrium geometry:
  real, intent(out), dimension(ntheta,nsmax) :: Rarr   ! R(theta,rho)
  real, intent(out), dimension(ntheta,nsmax) :: Zarr   ! Z(theta,rho)
!
!  normalized sqrt(toroidal flux):
  real, intent(out), dimension(nsmax) :: rho           ! sqrt(phi)/sqrt(philim)
!
!  poloidal flux profile
  real, intent(out), dimension(nsmax) :: psi           ! Webers/radian
!
!  MHD pressure
  real, intent(out), dimension(nsmax) :: pmhd          ! (see units_option)
!
!  q profil
  real, intent(out), dimension(nsmax) :: qmhd
!
!  g profile
  real, intent(out), dimension(nsmax) :: gmhd          ! (see units_option)
!
!  total enclosed flux
  real, intent(out) :: tflux                           ! Webers
!
!  total toroidal plasma current
  real, intent(out) :: pcur                            ! Amps
!
!----> completion code
!
  integer, intent(out) :: ierr     ! completion code, 0=OK
!
!----------------------------------------------------------------------
!
!  local variables and arrays
!
  real :: zconvrz          ! R,Z units conversion factor
  real :: zconvg           ! R*Bt conversion factor
!
  integer :: imax          ! nsmax-1
  real, dimension(nsmax-1) :: zbuf      ! profile buffer
  integer :: istype        ! profile x axis type
!
  integer :: isym          ! updown (a)symmetry flag
!
  real :: bzxr             ! R*Bt (vacuum)
  real :: q0               ! q on axis
!
  character(10) :: ztest   ! test "units_option"
!
  character(64) :: zlbl
  character(32) :: zuns
!
  real :: raxis,zaxis      ! mag. axis location
  integer :: iraxis,izaxis ! flag if raxis/zaxis scalars unavailable
!
  integer :: iadr          ! function index
!
  integer igot,istat,inx,inxp1
  integer lunzer
!----------------------------------------------------------------------
!
  ztest = units_option
  call uupper(ztest)
  if(ztest.eq.'TRANSP') then
     zconvrz=1.0
     zconvg=1.0
  else if(ztest.eq.'MKS') then
     zconvrz=0.01
     zconvg=0.01
  else
     call zermsg(' ?t1mhdeq:  invalid value of "units_option":  '// &
          units_option)
     go to 900
  endif
!
!-----------------------
!  read scalars
!
  write(lunzer(0),*) ' %t1mhdeq:  reading scalars...'
!
  call t1scalar('BZXR',zlbl,zuns,ztime,zdelta, bzxr, ierr)
  if(ierr.ne.0) then
     call zermsg(' ?t1mhdeq:  "BZXR" read failed.')
     go to 900
  endif
  bzxr=zconvg*bzxr   ! used to form g profile (R*Bt profile).
!
  call t1scalar('PCUR',zlbl,zuns,ztime,zdelta, pcur, ierr)
  if(ierr.ne.0) then
     call zermsg(' ?t1mhdeq:  "PCUR" read failed.')
     go to 900
  endif
!
  call t1scalar('TFLUX',zlbl,zuns,ztime,zdelta, tflux, ierr)
  if(ierr.ne.0) then
     call zermsg(' ?t1mhdeq:  "TFLUX" read failed.')
     go to 900
  endif
!
  call t1scalar('Q0',zlbl,zuns,ztime,zdelta, q0, ierr)
  if(ierr.ne.0) then
     call zermsg(' ?t1mhdeq:  "Q0" read failed.')
     go to 900
  endif
!
  iraxis=0
  call t1scalar('RAXIS',zlbl,zuns,ztime,zdelta, raxis, ierr)
  if(ierr.ne.0) then
     call zermsg(' %t1mhdeq:  "RAXIS" read failed, fixup applied.')
     iraxis=1
     ierr=0
  endif
!
  izaxis=0
  call t1scalar('YAXIS',zlbl,zuns,ztime,zdelta, zaxis, ierr)
  if(ierr.ne.0) then
     call zermsg(' %t1mhdeq:  "YAXIS" read failed, fixup applied.')
     izaxis=1
     ierr=0
  endif
!
!-----------------------
!  read profiles
!
  write(lunzer(0),*) ' %t1mhdeq:  reading profiles...'
!
  imax=nsmax-1
!
  call t1profil('XB',zlbl,zuns,ztime,zdelta, istype, zbuf, imax, igot, ierr)
  if(ierr.ne.0) then
     call zermsg(' ?t1mhdeq:  "XB" profile read failed.')
     go to 900
  endif
  inx=igot
  inxp1=inx+1
  rho(1)=0.0
  rho(2:inxp1)=zbuf(1:inx)
!
  nsgot = inxp1         ! **** number of surfaces, including axis.
!
  call t1profil('PLFLX',zlbl,zuns,ztime,zdelta, istype, zbuf, imax, igot, ierr)
  if(ierr.ne.0) then
     call zermsg(' ?t1mhdeq:  "PLFLX" profile read failed.')
     go to 900
  endif
  if(igot.ne.inx) then
     write(lunzer(0),*) ' ?t1mhdeq:  "PLFLX" & "XB" size inconsistency', &
          ' inx = ',inx, ' igot = ',igot
     go to 900
  endif
  psi(1)=0.0
  psi(2:inxp1)=zbuf(1:inx)
!
  call t1profil('PMHD_IN',zlbl,zuns,ztime,zdelta, istype, zbuf, imax, igot, &
       ierr)
  if(ierr.ne.0) then
     call zermsg(' %t1mhdeq:  "PMHD_IN" not found, reverting to "PTOWB".')
     call t1profil('PTOWB',zlbl,zuns,ztime,zdelta, istype, zbuf, imax, igot, &
          ierr)
     if(ierr.ne.0) then
        call zermsg(' ?t1mhdeq:  "PTOWB" profile read failed.')
        go to 900
     endif
  endif
  if(igot.ne.inx) then
     write(lunzer(0),*) &
          ' ?t1mhdeq:  "PMHD_IN or PTOWB" & "XB" size inconsistency', &
          ' inx = ',inx, ' igot = ',igot
     go to 900
  endif
!
!  shift to bdy's
!
  pmhd(1)=zbuf(1)
  pmhd(2:inx)=0.5*(zbuf(1:inx-1)+zbuf(2:inx))
  pmhd(inxp1)=max(0.1*zbuf(inx),(zbuf(inx)-(pmhd(inx)-zbuf(inx))))
!
  call t1profil('Q',zlbl,zuns,ztime,zdelta, istype, zbuf, imax, igot, ierr)
  if(ierr.ne.0) then
     call zermsg(' ?t1mhdeq:  "Q" profile read failed.')
     go to 900
  endif
  if(igot.ne.inx) then
     write(lunzer(0),*) ' ?t1mhdeq:  "Q" & "XB" size inconsistency', &
          ' inx = ',inx, ' igot = ',igot
     go to 900
  endif
  qmhd(1)=q0
  qmhd(2:inxp1)=zbuf(1:inx)
!
  call t1profil('GFUN',zlbl,zuns,ztime,zdelta, istype, zbuf, imax, igot, ierr)
  if(ierr.ne.0) then
     call zermsg(' ?t1mhdeq:  "GFUN" profile read failed.')
     go to 900
  endif
  if(igot.ne.inx) then
     write(lunzer(0),*) ' ?t1mhdeq:  "GFUN" & "XB" size inconsistency', &
          ' inx = ',inx, ' igot = ',igot
     go to 900
  endif
!
  gmhd(2:inxp1)=zbuf(1:inx)
  gmhd(1)=gmhd(2)+0.33333333*(gmhd(2)-gmhd(3))
  gmhd=gmhd*bzxr        ! convert to R*Bt (units already converted)
!
!-----------------------
! read equilibrium...
!
  call rp_geteq(ztime,zdelta,inx,ntheta,theta,Rarr(1,2),Zarr(1,2),ierr)
  if(ierr.ne.0) go to 900
!
  if(iraxis.eq.0) then
     Rarr(1:ntheta,1)=raxis
  else
     call t1mhdeq_axtrap(Rarr,ntheta)
  endif
!
  if(izaxis.eq.0) then
     Zarr(1:ntheta,1)=zaxis
  else
     call rp_eq_symflag(isym)
     if(isym.eq.1) then
        call t1mhdeq_axtrap(Zarr,ntheta)
     else
        Zarr(1:ntheta,1)=0.0
     endif
  endif
!
  Rarr = zconvrz * Rarr
  Zarr = zconvrz * Zarr
!
  return                      ! normal exit
!----------------------------------------------------------------------
!
900 continue
  Rarr = 0.0
  Zarr = 0.0
  rho = 0.0
  psi = 0.0
  pmhd = 0.0
  qmhd = 0.0
  gmhd = 0.0
  tflux = 0.0
  pcur = 0.0
!
  ierr = 1
  return                      ! error exit
!
end subroutine t1mhdeq
 
subroutine t1mhdeq_axtrap(rzarr,ntheta)
!
  implicit none
!
!  extrapolate axis from 1st 2 surfaces
!
  integer, intent(in) :: ntheta                      ! no. of theta pts / surf.
  real, intent(inout), dimension(ntheta,3) :: rzarr  ! 1st 3 surfaces...
!
! rzarr can be R or Z.  rzarr(1:ntheta,2:3) are known on input
!                       rzarr(1:ntheta,1) is computed on output
!
!------------------------------
  real rzavg2,rzavg3
!------------------------------
!
  rzavg2=sum(rzarr(1:ntheta,2))/ntheta
  rzavg3=sum(rzarr(1:ntheta,3))/ntheta
!
  rzarr(1:ntheta,1)=(4.0*rzavg2 - rzavg3)/3.0
!
  return
  end
