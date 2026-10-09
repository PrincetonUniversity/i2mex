subroutine r8bsmoo(zprof,zoption)
  !
  !  smooth radial profile -- triangular smooth
  !  (sum of cubics i.e. definite integrals of products of linear pieces).
  !
  !     *** smoothing half-width (in xi):  namelist control DXBSMOO ***
  !
  !  weighting fcn about each point x0 is:
  !
  !   w(x)=
  !     zero for |x-x0| >= dxbsmoo
  !     (1/dxbsmoo)-|x-x0|/(dxbsmoo*dxbsmoo) for |x-x0| <= dxbsmoo
  !
  !  this makes a symmetrical triangle of area 1 centered on x0 and with base
  !  2*dxbsmoo wide, height 1/dxbsmoo.
  !
  !  the smoothed function is treated as a sequence of linear segments.
  !
  !  f(x) = f(j)-x(j)*S + x*S    ... for x btw x(j) and x(j+1)
  !                                  and S = (f(j+1)-f(j))/(x(j+1)-x(j))
  !
  !  piecewise integration of w(x)*f(x) (quadratic pieces) btw data pts
  !  yields a sum of cubics evaluated as definite integrals.
  !
  !  curve will be "normalized" via znorm factor, prior to smooth, then
  !  "un-normalized" on the way back out.  Allows user to control the "units"
  !  of the space in which the smooth operation is carried out
  !
  !  zprof(...) -- profile to be smoothed
  !  zoption(...) -- if zoption(lcentr).gt.0.0 then
  !                    ...copy whole array into znorm
  !                    ...use dxbsmoo as smoothing parameter
  !                  if zoption(lcentr).le.0.0 then
  !                    ...set znorm to 1.0
  !                    ...set iopt = -zoption(lcentr)
  !                       ...iopt = 0 means standard bdy conditions
  !                       ...iopt = 1 means:  use reflection to fix point
  !                          at outer bdy, position lep1.
  !                    ...use zoption(lcp1) as smoothing parameter
  !
  !
  !
  !----------
  !
  use iso_c_binding, only: fp => c_double
  use r8bsmoo_mod
  !
  implicit none
  integer :: iopt,j,inzp1,in3zp1,iz
  integer :: ii,icen,ism,iamdone
  !============
  ! idecl:  explicitize implicit REAL declarations:
  real(fp) :: zxsmoo,zedge,zdpt,zdiff,zdbxi,zdbx2,zxcen
  real(fp) :: zx1b,zx2b,zf1b,zf2b,zx1a,zx2a,zf1a,zf2a,zfac,zdelx
  real(fp) :: zsumx,zsmx2,zslope,zcns,za,zb,zc,zans
  !============
  real(fp) :: zprof(MJ)                    ! in:  unsmoothed, out:  smoothed
  real(fp) :: zoption(MJ)                  ! controls
  !
  real(fp) :: znorm(MJ)                    ! presmooth normalization
  !
  real(fp) :: zwork(MJ)                    ! workspace
  real(fp) :: zworkx(3*MJ)
  real(fp) :: zworkf(3*MJ)
  !
  !-----------------------------
  ! check the smoothing parameter
  !
  if(zoption(lcentr).gt.0.0_fp) then
    ! beam code output smooth
    if(dxbsmoo.le.0.0_fp) return
    dxbsmoo=max(0.01_fp,min(0.2_fp,dxbsmoo))
    zxsmoo=dxbsmoo
    iopt=0
    do j=1,mj
      znorm(j)=zoption(j)
    end do
  else
    ! NCLASS code input smooth
    iopt=-zoption(lcentr)
    zxsmoo=zoption(lcp1)
    if(zxsmoo.le.0.0_fp) return
    zxsmoo=max(0.01_fp,min(0.2_fp,zxsmoo))
    do j=1,mj
      znorm(j)=1.0_fp
    end do
  end if
  !
  !----------------
  !
  !  boundary conditions:
  !    extrapolation by symmetric reflection at center
  !    flat extrapolation at edge
  !
  inzp1=nzones+1
  in3zp1=3*nzones+1
  do iz=1,nzones
    j=iz+lcentr-1
    zworkf(iz+nzones)=zprof(j)/znorm(j)
    zworkf(inzp1-iz) =zprof(j)/znorm(j)     ! reflection
    zworkf(in3zp1-iz)=zprof(ledge)/znorm(ledge)  ! flat extrap.
    zworkx(iz+nzones)=xi(j,2)
    zworkx(inzp1-iz)=-xi(j,2)
    zworkx(in3zp1-iz)=2*xi(lep1,1)-xi(j,2)
  end do
  !
  if(iopt.eq.1) then
    zedge=zprof(lep1)/znorm(lep1)
    zworkf(2*nzones+1)=zedge
    !  anti-symmetric reflection beyond edge point; this causes the
    !  symmetric smoothing operator to fix the edge point.
    do ii=2,nzones
      j=ledge+2-ii
      zdpt=zprof(j)/znorm(j)
      zdiff=zdpt-zedge
      zworkf(2*nzones+ii)=zedge-zdiff
    end do
  end if
  !
  !  working from zworkf & zworkx compute smoothed results and store
  !  back in zprof
  !
  zdbxi=1.0_fp/zxsmoo
  zdbx2=zdbxi*zdbxi
  !
  do iz=1,nzones
    j=iz+lcentr-1
    zprof(j)=0.0_fp  ! clear the sum
    icen=iz+nzones
    zxcen=zworkx(icen)
    zx1b=0.0_fp
    zx2b=0.0_fp
    zf1b=zworkf(icen)
    zf2b=zf1b
    do ism=1,nzones
      zx1a=zx1b
      zx2a=zx2b
      zf1a=zf1b
      zf2a=zf2b
      zx1b=zxcen-zworkx(icen-ism) ! reflected, for convenience
      zx2b=zworkx(icen+ism)-zxcen
      zf1b=zworkf(icen-ism)
      zf2b=zworkf(icen+ism)
      iamdone=0
      !
      !  deal with crossing left edge of weighting triangle
      if(zx1b.ge.zxsmoo) then
        iamdone=iamdone+1
        zfac=(zxsmoo-zx1a)/(zx1b-zx1a)
        zf1b=zf1a+(zf1b-zf1a)*zfac
        zx1b=zxsmoo
      end if
      !
      !  deal with crossing right edge of weighting triangle
      if(zx2b.ge.zxsmoo) then
        iamdone=iamdone+1
        zfac=(zxsmoo-zx2a)/(zx2b-zx2a)
        zf2b=zf2a+(zf2b-zf2a)*zfac
        zx2b=zxsmoo
      end if
      !
      !  add in left piece
      if(zx1a.lt.zxsmoo) then
        zdelx=zx1b-zx1a
        zsumx=zx1b+zx1a
        zsmx2=zx1b*zx1b+zx1b*zx1a+zx1a*zx1a
        zslope=(zf1b-zf1a)/zdelx
        zcns=(zf1a-zx1a*zslope)
        za=zcns*zdbxi
        zb=0.5_fp*(zslope*zdbxi-zcns*zdbx2)
        zc=-0.33333333333333333_fp*zslope*zdbx2
        zans=zdelx*(za+ zsumx*zb +zsmx2*zc)
        zprof(j)=zprof(j)+zans
      end if
      !
      !  add in right piece
      if(zx2a.lt.zxsmoo) then
        zdelx=zx2b-zx2a
        zsumx=zx2b+zx2a
        zsmx2=zx2b*zx2b+zx2b*zx2a+zx2a*zx2a
        zslope=(zf2b-zf2a)/zdelx
        zcns=(zf2a-zx2a*zslope)
        za=zcns*zdbxi
        zb=0.5_fp*(zslope*zdbxi-zcns*zdbx2)
        zc=-0.33333333333333333_fp*zslope*zdbx2
        zans=zdelx*(za+ zsumx*zb +zsmx2*zc)
        zprof(j)=zprof(j)+zans
      end if
      !
      !  test if done with this triangle
      if(iamdone.eq.2) go to 100
    end do                          ! ism, triangle piece loop
100 continue
    zprof(j)=zprof(j)*znorm(j)     ! un-normalize
  end do                             ! iz, profile zone loop
  !
  !  all done
  !
  return
end subroutine r8bsmoo
