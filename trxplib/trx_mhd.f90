subroutine trx_mhd(ntheta,ierr)
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
  integer i,iR,iZ,igot,idum,idumR,idumZ
  integer :: lunzer
!
  real Rarr(ntheta,nsurf),Zarr(ntheta,nsurf)
  real*8 r8Rarr(ntheta,nsurf),r8Zarr(ntheta,nsurf)
!
  real zrho(nsurf),zpsi(nsurf),zpmhd(nsurf),zqmhd(nsurf),zgmhd(nsurf)
  real*8 zbuf(nsurf)
  real ztflux,zpcur
!
  real*8, dimension(:), allocatable :: arr1
  real*8 :: zrcen1,zrcen2,zzcen1,zzcen2,zrtest,zztest,testval,testval2,zdpsi
  real*8 :: zmina(3),zresult
  integer :: ircen1,ircen2,izcen1,izcen2,iraxis,izaxis
  logical :: iexist
!----------------------------
!
  call trx_ready('trx_mhd',ierr)
  if(ierr.ne.0) return
!
  zpi=acos(-1.0E0_R8)
!
  do i=1,ntheta
     ztheta(i) = -zpi + (i-1)*2.0E0_R8*zpi/(ntheta-1)
     if(th_reverse.eq.0) then
        r4ztheta(i)=ztheta(i)
     else
        r4ztheta(ntheta-i+1)=ztheta(i)
     endif
  enddo
!
  call eqm_chi(ztheta,0,ntheta,1.0E-4_R8,id_chi,ierr)
  if(ierr.ne.0) return
!
  call t1mhdeq(time0,delta_t,nsurf,ntheta,igot,'MKS',r4ztheta, &
       Rarr,Zarr,zrho,zpsi,zpmhd,zqmhd,zgmhd,ztflux,zpcur,ierr)
  if(ierr.ne.0) return
!
  ztime0=time0   ! R8
  zdelt =delta_t ! R8

  ! make Psi profiles zero on axis -- they should be close already

  if(abs(zpsi(1)).gt.1.0d-4*abs(zpsi(nsurf))) then
     write(lunzer(0),*) ' %trx_mhd: 1d profile Psi(rho): Psi(0)=',zpsi(1)
     write(lunzer(0),*) '  reset to zero...'
  endif
  zdpsi = zpsi(1)
  zpsi = zpsi - zdpsi
!
  if((nRfree.gt.0).and.(nZfree.gt.0)) then
     ! retrieve Psi(R,Z) also...
     allocate(arr1(nRfree*nZfree))

     call rpexist_profile('PSIRZ',iexist)

     if (iexist) then
        call r8_t1profil('PSIRZ',zlbl,zuns,ztime0,zdelt, &
             idum,arr1,nRfree*nZfree,igot,ierr)

     else if (ifound_psi0.eq.1) then
        call r8_t1scalar('PSI0_TR',zlbl,zuns,ztime0,zdelt,Psi0_mhd,ierr)

        if (ierr==0) then
           call r8_t1profil('APSIRZ',zlbl,zuns,ztime0,zdelt, &
                idum,arr1,nRfree*nZfree,igot,ierr)
           if (ierr==0) then
              arr1 = arr1 - Psi0_mhd  ! psi relative to machine axis now psi relative to mag axis
           end if
        end if
     else
        ierr=1
     end if

     if(ierr.ne.0) then
        deallocate(arr1)
        return
     endif
     iZ=1
     iR=0
     do i=1,nRfree*nZfree
        iR=iR+1
        if(iR.gt.nRfree) then
           iR=1
           iZ=iZ+1
        endif
        psiRZ_free(iR,iZ)=arr1(i)
     enddo
     ! also need axis location -- cannot ignore error
     call r8_t1scalar('RAXIS',zlbl,zuns,ztime0,zdelt,raxis_mhd,ierr)
     if(ierr.ne.0) return
     raxis_mhd = 0.01_R8*raxis_mhd

     call r8_t1scalar('YAXIS',zlbl,zuns,ztime0,zdelt,zaxis_mhd,ierr)
     if(ierr.ne.0) return
     zaxis_mhd = 0.01_R8*zaxis_mhd

     if(ifound_psi0.eq.1) then
        call r8_t1scalar('PSI0_TR',zlbl,zuns,ztime0,zdelt,Psi0_mhd,ierr)
        if(ierr.ne.0) return
        call eqm_save_psi0(Psi0_mhd,ierr)
        if(ierr.ne.0) return
     else
        Psi0_mhd=0.0_R8
     endif

     ! use axis information to normalize Psi(R,Z): 0 on axis (approximately)
     zrcen1 = 0.9*raxis_mhd
     zrcen2 = 1.1*raxis_mhd
     zzcen1 = zaxis_mhd - 0.1*raxis_mhd
     zzcen2 = zaxis_mhd + 0.1*raxis_mhd

     zrtest = 100*raxis_mhd
     zztest = zrtest

     ircen1=0
     ircen2=0
     do ir=1,nRfree
        if(Rgrid_free(ir).lt.zrcen1) ircen1=ir
        if(Rgrid_free(ir).gt.zrcen2) then
           ircen2=ir
           exit
        endif
        testval = abs(Rgrid_free(ir)-Raxis_mhd)
        if(testval.lt.zRtest) then
           iraxis=ir
           zRtest=testval
        endif
     enddo

     izcen1=0
     izcen2=0
     do iz=1,nZfree
        if(Zgrid_free(iz).lt.zZcen1) izcen1=iz
        if(Zgrid_free(iz).gt.zZcen2) then
           izcen2=iz
           exit
        endif
        testval = abs(Zgrid_free(iz)-Zaxis_mhd)
        if(testval.lt.zZtest) then
           izaxis=iz
           zZtest=testval
        endif
     enddo

     ! see if Psi(axis) local minimum or maximum

     testval = PsiRZ_free(iraxis,izaxis)
     testval2 = (PsiRZ_free(ircen1,izcen1)+PsiRZ_free(ircen2,izcen1)+ &
          PsiRZ_free(ircen1,izcen2)+PsiRZ_free(ircen2,izcen2))*0.25_R8

     if(testval.gt.testval2) then
        ! local maximum: flip Psi

        PsiRZ_free = -PsiRZ_free
        testval = -testval
     endif

     ! find actual minimum
     do iZ=izcen1,izcen2
        do iR=ircen1,ircen2
           if(PsiRZ_free(iR,iZ).lt.testval) then
              testval = PsiRZ_free(iR,iZ)
              iRaxis=iR
              iZaxis=iZ
           endif
        enddo
     enddo

     ! 2nd order adjustment -- use CONTAIN'ed subroutine

     call fmin1d(Rgrid_free(iRaxis-1:iRaxis+1), &
          PsiRz_free(iRaxis-1:iRaxis+1,iZaxis-1), zmina(1))

     call fmin1d(Rgrid_free(iRaxis-1:iRaxis+1), &
          PsiRz_free(iRaxis-1:iRaxis+1,iZaxis), zmina(2))

     call fmin1d(Rgrid_free(iRaxis-1:iRaxis+1), &
          PsiRz_free(iRaxis-1:iRaxis+1,iZaxis+1), zmina(3))

     call fmin1d(Zgrid_free(iZaxis-1:iZaxis+1), zmina, zresult)

     testval = zresult

     if(abs(testval).gt.1.0d-4*abs(zpsi(nsurf))) then
        write(lunzer(0),*) ' %trx_mhd: 2d profile Psi(R,Z): Psi(axis)=', &
             testval
        write(lunzer(0),*) '  reset to zero...'
     endif

     PsiRZ_free = PsiRZ_free - testval
        
  else
     ! no Psi(R,Z): try to read axis values but ignore error
     call r8_t1scalar('RAXIS',zlbl,zuns,ztime0,zdelt,raxis_mhd,ierr)
     raxis_mhd = 0.01_R8*raxis_mhd
     call r8_t1scalar('YAXIS',zlbl,zuns,ztime0,zdelt,zaxis_mhd,ierr)
     zaxis_mhd = 0.01_R8*zaxis_mhd
     ierr = 0

  endif
!
  call rp_eq_symflag(ksym)  ! updown (a)symmetry
  call rp_eq_nmoms(nmoms)   ! no. of moments in TRANSP equilibrium
  call xmoments_kmom_set(nmoms)  ! tell XPLASMA
!
! set up 1d profiles
!
  zbuf=ztflux*zrho*zrho
  call eqm_rhofun(2,id_rho,'phitor',zbuf,1,0.0E0_R8,0,0.0E0_R8,idum,ierr)
  call trx_ustore(ierr,idum,id_rho,-99,'Wb',0)
  if(ierr.ne.0) return
!
  zbuf=zpsi
  call eqm_rhofun(2,id_rho,'psi',zbuf,1,0.0E0_R8,0,0.0E0_R8,idum,ierr)
  call trx_ustore(ierr,idum,id_rho,-99,'Wb/rad',0)
  if(ierr.ne.0) return
!
  zbuf=zgmhd
  call eqm_rhofun(2,id_rho,'g',zbuf,1,0.0E0_R8,1,0.0E0_R8,idum,ierr)
  call trx_ustore(ierr,idum,id_rho,-99,'T*m',0)
  if(ierr.ne.0) return
!
  zbuf=zpmhd
  call eqm_rhofun(2,id_rho,'pmhd',zbuf,1,0.0E0_R8,0,0.0E0_R8,idum,ierr)
  call trx_ustore(ierr,idum,id_rho,-99,'Pa',0)
  if(ierr.ne.0) return

  call eqm_mark_pmhd(idum,ierr)  ! mark the pressure profile
  if(ierr.ne.0) return
!
  zbuf=zqmhd
  call eqm_rhofun(2,id_rho,'q',zbuf,0,0.0E0_R8,0,0.0E0_R8,idum,ierr)
  call trx_ustore(ierr,idum,id_rho,-99,' ',0)
  if(ierr.ne.0) return
!
  r8Rarr=Rarr
  r8Zarr=Zarr
!
  call eqm_rzmag(r8Rarr,r8Zarr,ntheta,nsurf,2*isign_rzbc,idumR,idumZ,ierr)
  call trx_ustore(ierr,idumR,0,-99,'m',0)
  call trx_ustore(ierr,idumZ,0,-99,'m',0)
  if(ierr.ne.0) return
!
  return

CONTAINS
  subroutine fmin1d(xa,fa,ans)

    real*8, dimension(:) :: xa   ! x array -- 3 values actually
    real*8, dimension(:) :: fa   ! f(x) array -- 3 values

    real*8, intent(out) :: ans   ! min value

    ! use a parabolic fit to 3 pts to estimate min/max value
    ! (in trx_mhd it will be a minimum, because the search is in vicinity
    ! of minimum of Psi(R,Z)).

    !-------------------------------------------
    real*8 :: xm,xp   ! distances from x(2)
    real*8 :: fm,fp   ! f(1)-f(2),f(3)-f(2)
    real*8 :: a,b,denom,xminloc
    !-------------------------------------------
    ! check against singular cases...

    if((fa(1).eq.fa(2)).AND.(fa(2).eq.fa(3))) then
       ans = fa(2)
       return
    endif

    if((fa(1).lt.fa(2)).AND.(fa(2).lt.fa(3))) then
       ans = fa(1)
       return
    endif

    if((fa(1).gt.fa(2)).AND.(fa(2).gt.fa(3))) then
       ans = fa(3)
       return
    endif

    !-------------------------------------------
    ! OK do parabolic fit

    xm = xa(2)-xa(1)
    xp = xa(3)-xa(2)

    fm = fa(1)-fa(2)
    fp = fa(3)-fa(2)

    ! look for parabola p which satisfies: p(0)=0, p(-xm)=fm, p(xp)=fp
    !   p = a*x**2 + b*x;  a*xm**2 - b*xm = fm;  a*xp**2 + b*x = fp

    denom = xm*xp*(xm+xp)
    a = (fm*xp + fp*xm)/denom
    b = (fp*xm*xm - fm*xp*xp)/denom

    !   find p'=0,  p' = 2*a*x + b = 0 => x_min = -b/2*a;
    !        then p_min = a*x_min**2 + b*x_min

    xminloc = -b/(2*a)

    ans = xminloc*(a*xminloc + b)

    ans = ans + fa(2)   ! add this back in...

  end subroutine fmin1d

end subroutine trx_mhd
 
!----------------------------------------------------------------------
!  the following is used for testing...
 
subroutine trx_set_threverse(ival)
 
  use trx_module
  implicit NONE
 
  integer ival            ! =1:  reverse theta; =0: don't.
 
  th_reverse=max(0,min(1,ival))
 
  return
 
end subroutine trx_set_threverse
