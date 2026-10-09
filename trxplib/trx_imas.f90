subroutine trx_imas(ntheta,ierr)
  use trx_module
  implicit NONE
!
!  read TRANSP MHD equilibrium data; store in xplasma for use via
!  interpolation routines
!
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
  integer, intent(IN) :: ntheta        ! no. of theta grid pts [-pi,pi]
  integer, intent(OUT) :: ierr         ! completion code, 0=OK
!
!----------------------------
!
  real*4 r4ztheta(ntheta)
  real*8 ztheta(ntheta),zpi
!
! for Psi(R,Z) acquisition (free bdy mode)
  character*64 zlbl
  character*32 zuns
  real*8 :: ztime0,zdelt
!
  integer igot,idum
  integer :: lunzer
!
  real Rarr(ntheta,nsurf),Zarr(ntheta,nsurf)
  real*8 r8Rarr(ntheta,nsurf),r8Zarr(ntheta,nsurf)
!
  real zrho(nsurf),zpsi(nsurf),zpmhd(nsurf),zqmhd(nsurf),zgmhd(nsurf)
  real*8 zbuf(nsurf)
  real ztflux,zpcur
  integer isurf,inx,inxp1
!
  real*8, dimension(:), allocatable :: arr1
  real*8 :: zrcen1,zrcen2,zzcen1,zzcen2,zrtest,zztest,testval,testval2,zdpsi
  real*8 :: zmina(3),zresult
  integer :: ircen1,ircen2,izcen1,izcen2,iraxis,izaxis
  logical :: iexist
!----------------------------
!
  call trx_ready('trx_imas',ierr)
  if(ierr.ne.0) return
!
  ztime0=time0   ! R8
  zdelt =delta_t ! R8
  isurf=nsurf-1
  allocate(arr1(isurf))

  call rpexist_profile('CUR',iexist)

  if (iexist) then
     call r8_t1profil('CUR',zlbl,zuns,ztime0,zdelt, &
          idum,arr1,isurf,igot,ierr)
     if(ierr.ne.0) then
        call zermsg(' ?trx_imas:  "CUR" profile read failed.')
        ierr=1
     endif
  else
     call zermsg(' ?trx_imas:  "CUR" profile does not exist.')
     ierr=1
  end if
  if(ierr.ne.0) then
     deallocate(arr1)
     return
  endif

  inx=igot
  inxp1=inx+1

!
! set up 1d profiles
!
  zbuf=0.0
  zbuf(2:inxp1)=arr1(1:inx)*1.0e-4 ! A/cm**2 -> A/m**2
  call eqm_rhofun(2,id_rho,'icurt',zbuf,1,0.0E0_R8,0,0.0E0_R8,idum,ierr)
  call trx_ustore(ierr,idum,id_rho,-99,'A/m**2',0)
  if(ierr.ne.0) then
     deallocate(arr1)
     return
  endif
!
  call rpexist_profile('PLJB',iexist)
  arr1=0.0
  if (iexist) then
     call r8_t1profil('PLJB',zlbl,zuns,ztime0,zdelt, &
          idum,arr1,isurf,igot,ierr)
     if(ierr.ne.0) then
        call zermsg(' ?trx_imas:  "PLJB" profile read failed.')
        ierr=1
     endif
  else
     call zermsg(' ?trx_imas:  "PLJB" profile does not exist.')
     ierr=1
  end if
  if(ierr.ne.0) then
     deallocate(arr1)
     return
  endif
!
  zbuf=0.0
  zbuf(2:inxp1)=arr1(1:inx)*1.0e-4 ! A*T/cm**2 -> A*T/m**2
  call eqm_rhofun(2,id_rho,'curpll',zbuf,1,0.0E0_R8,0,0.0E0_R8,idum,ierr)
  call trx_ustore(ierr,idum,id_rho,-99,'A*T/m**2',0)
  if(ierr.ne.0) then
     deallocate(arr1)
     return
  endif
!
  return
end subroutine trx_imas
 
