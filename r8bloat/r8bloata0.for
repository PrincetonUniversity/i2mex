      subroutine r8bloata0(symflag,lcentr,nzones,xiblo,inbblo,mtbl,mj,
     >                  rmc,ymc,nmom,mimom,nonlin,
     >                  rmcx,ymcx,ierr)
C
C  bloat subroutine -- updown *asymmetric* -- calling arguments, no COMMON!
C  dmc 27 Dec 1996
C
C  algorithm rewritten -- DMC Aug 14 2003 -- continuing reliability issues
C  with old algorithm not to be overcome.
C
C  Restriction:  thetaA(theta) = atan2((R(theta)-Rcen),(Z(theta)-Zcen))
C  must be monotonic with some margin: d(thetaA)/d(theta) > 0.01D0 ...
C
C  In other words, any straight line drawn out from the magnetic axis
C  must only cross the boundary surface once, and also cannot be
C  tangent to the boudary surface.
C
C  this routine copies the inside-the-plasma moments set to the bloated
C  moments set, and uses a numerical algorithm to create well behaved
C  surfaces beyond the plasma bdy.  these surfaces provide a coordinate
C  system with reasonable properties (|J|>0 always) and a Fourier moments
C  representation, but the ones beyond the boundary are not real magnetic
C  surfaces!
C
C  The (rho,theta)->(R,Z) mapping is continuous but not differentiable
C  accross the plasma boundary; separate splines should be used for the
C  core plasma and vacuum regions.
C
C input:
C  symflag -- TRUE if input equilibrium is updown symmetric (even though
C             the moments arrays have slots for assymetric moment terms).
C
C  lcentr -- index to centermost surface (mag. axis point).
C            typically 1; 2 (TRCOM's LCENTR) when called from TRANSP.
C  nzones -- the number of flux zones in the plasma (the no. of surfaces
C            counting the mag. axis is nzones+1).
C
C  xiblo(..)  -- x grid for extrapolated space beyond plasma bdy
C  inbblo -- no. of points in use in extrapolated grid.
C  mtbl   -- max no. of grid pts (array dimension for xiblo, etc...)
C  mj     -- max no. of grid pts *inside plasma* (array dim for rmom0, etc.)
C
C  rmc(...) -- inside the plasma R moments
C  ymc(...) -- inside the plasma Y moments
C
C   caution:  dimension ordering follows TRANSP COMMON!
C
C  nmom   -- no. of Fourier moments in use
C  mimom  -- max no. of Fourier moments (array dimension)
C  nonlin -- LUN for writing error message
C
C output:
C
C   caution:  dimension ordering follows TRANSP COMMON!
C
C  rmcx(...) -- R moments (extrapolated grid)
C  ymcx(...) -- Y moments (extrapolated grid)
C
C  ierr -- status code: 0 = OK
C
C------------------------------------------------------------------
C
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
C
C  passed items:
      INTEGER nzones,inbblo,mtbl,mj,nmom,mimom,nonlin,lcentr
C
      logical symflag
C
      REAL*8 xiblo(mtbl)
C
      REAL*8 rmc(mj,0:mimom,2),ymc(mj,0:mimom,2)
C
      REAL*8 rmcx(mtbl,0:mimom,2)
      REAL*8 ymcx(mtbl,0:mimom,2)
C
      integer ierr
C
C------------------------------------------------------------------
c
      include 'r8tkbparm.inc'
c
C------------------------------------------------------------------
C  local items:
C
      integer inzp1,inzp2,lcm1
      integer j,jzb,im,ics,jbx0,jbxm1,itk,jnext,jmax
      real*8 zshift
C
C  Axis:
      real*8 Rcen,Zcen
C  ThetaA(theta) -- monotonic increasing & invariant
      real*8 Theta(ntk+1),ThetaA(0:ntk+1),rcoef(ntk),zcoef(ntk)
C
C  axis to surface distance -- current surface & previous surface
C
      real*8 dist(ntk),distp(ntk),d_dist(ntk)
C
C  surface contour
      real*8 Rpts(ntk),Zpts(ntk)
C
C  average distance; average change in distance from prev. surface
      real*8 distav,d_distav
C
C  locally smoothed distance
C
      real*8 distsm(ntk)
C
      real*8, parameter :: c2pi = 6.2831853071795862D+00
      real*8, parameter :: cpi =  3.1415926535897931D+00
C
C  smoothing parameter (for distsm):
C
      real*8, parameter :: delsm = c2pi/25
C
C  current step, relative to bdy xiblo value
C
      real*8 delta_x
C
C  surface expansion controls
C
      real*8 alpha,beta,gamma,const
      real*8 afac(ntk),bfac(ntk),gfac(ntk),pfac(ntk)
C
      real*8, parameter :: ZERO = 0.0_R8
      real*8, parameter :: ONE = 1.0_R8
C
      real*8, parameter :: zdthmin = 0.0100_R8  ! min d(ThetaA)/d(Theta).
C
C  fft
C
      real*8 zs1d(ntk),zc1d(ntk)
C
      character*60 mnlab(3)
C
C------------------------------------------------------------------
C   surface expansion algorithm:
C     (A) theta lines extend straight out from the boundary, colinear
C         from a central point (Rcen,Zcen).
C     (B) each successive flux surface is extended from the prior
C         pair of surfaces such as to
C           0.  grow a set of well behaved nested surfaces out from
C               the plasma boundary. (const)
C           1.  reduce global variations in distance from (Rcen,Zcen)
C               towards an average value. (alpha)
C           2.  reduce global variations in spacing of successive
C               surfaces towards an average value. (beta)
C           3.  smooth out local variations in distance from (Rcen,Zcen).
C               (gamma)
C------------------------------------------------------------------
C   first copy the interior moments
C
C   NO. OF BDY'S INSIDE PLASMA; no. of zones PLUS 1
C
      INZP1=NZONES+1
      INZP2=NZONES+2
C
      LCM1=LCENTR-1
      jmax=inbblo+lcm1
C
      DO 60 J = 1, INZP1
         JZB=J+LCM1
         DO IM = 0, NMOM
            DO ICS = 1, 2
               RMCX(JZB,IM,ICS) = RMC(JZB,IM,ICS)
               YMCX(JZB,IM,ICS) = YMC(JZB,IM,ICS)
            enddo
         enddo
 60   CONTINUE
C
C  form & check theta arrays
C
C  boundary surface
      JBX0=INZP1+LCM1
C
      do itk=1,ntk
         Theta(itk)=(itk-1)*c2pi/ntk
      enddo
C
C  get a central point
C
      call get_RZcen(jbx0)
C
C  distance from point
C
      call get_RZdist(jbx0)
 
      zshift=0
      ierr=0
c
c  these items are computed just once at the plasma boundary surface...
c
      do itk=1,ntk
         ThetaA(itk)=atan2(Zpts(itk)-Zcen,Rpts(itk)-Rcen) + zshift
         Rcoef(itk)=(Rpts(itk)-Rcen)/dist(itk)
         Zcoef(itk)=(Zpts(itk)-Zcen)/dist(itk)
         if(itk.gt.1) then
            if(ThetaA(itk).lt.ThetaA(itk-1)-cpi) then
               ThetaA(itk)=ThetaA(itk) + c2pi
               zshift = zshift + c2pi
               if(zshift.gt.1.1_R8*c2pi) ierr=ierr+1
            endif
            if(ThetaA(itk).lt.ThetaA(itk-1)+zdthmin*c2pi/ntk) then
#ifdef __DEBUG
               write(6,*) itk-1,' dth=',thetaA(itk)-thetaA(itk-1), 
     >              ' min=',zdthmin*c2pi/ntk
#endif
               ierr=ierr+1
            endif
         endif
      enddo
      if(ThetaA(1).lt.ThetaA(ntk)-c2pi+zdthmin*c2pi/ntk) then
#ifdef __DEBUG
         write(6,*) ntk,' dth=',thetaA(1)-thetaA(ntk)-c2pi,
     >        ' min=',zdthmin*c2pi/ntk
#endif
         ierr=ierr+1
      endif
      ThetaA(0)=ThetaA(ntk)-c2pi
      ThetaA(ntk+1)=ThetaA(1)+c2pi
C
      if(ierr.ne.0) then
         ierr=99                        ! theta not monotonic
         return
      endif
C
      Theta(ntk+1)=Theta(1)+c2pi
C
cdbg      mnlab(1)='r8bloats0 debug plot'
cdbg      mnlab(2)=' '
cdbg      mnlab(3)='R(theta),Z(theta)'
cdbg      call r8_grafx2(theta,Rpts,Zpts,ntk,'radians','m or cm',
cdbg     >   mnlab(1),mnlab(2),mnlab(3))
C
cdbg      mnlab(3)='ArcTan[(Z(theta)-Zcen),(R(theta)-Rcen)] vs. theta'
cdbg      call r8_grafx2(theta,thetaA,theta,ntk,'radians','radians',
cdbg     >   mnlab(1),mnlab(2),mnlab(3))
C
C	DO THE EXTRAPOLATION INTO THE BLOATED ZONES
C
C  one inside from boundary surface
      JBXM1=JBX0-1
C  find intercepts of lines from axis to boundary with JBXM1 surface...
C
      call find_distm1(JBXM1)
C
C  for each surface...
C
      jnext=jbx0
C
ctmp      alpha=2.0_R8
ctmp      call ureadr8(' alpha (for <D>/D factor):$',alpha)
ctmp      alpha=max(0.0_r8,min(100.0_r8,alpha))
 
ctmp      beta=3.0_R8
ctmp      call ureadr8(' beta (for <D-Dprev>/(D-Dprev) factor):$',beta)
ctmp      beta=max(0.0_r8,min(100.0_r8,beta))
 
ctmp      gamma=2.0_R8
ctmp      call ureadr8(' gamma (for Dsm/D factor):$',gamma)
ctmp      gamma=max(0.0_r8,min(100.0_r8,gamma))
 
ctmp      const=1.0_R8
ctmp      call ureadr8(' growth constant:$',const)
ctmp      const=max(0.5_r8,min(2.0_r8,const))
C
      do
         jnext=jnext+1
         if(jnext.gt.jmax) exit
C
         delta_x=(xiblo(jnext)-xiblo(jnext-1))/
     >      (xiblo(jbx0)-xiblo(lcentr))
C
         d_dist = dist-distp
C
         call get_avdist
C
         call sm_dist
C
cdbg         mnlab(1)='r8bloats0 debug plot'; mnlab(2)=' '
cdbg         write(mnlab(2),1001) distav,d_distav
cdbg 1001    format('avg(dist)=',1pe11.4,', avg(d_dist)=',1pe11.4)
cdbg         mnlab(3)='dist,distp,distsm,d_dist'
cdbg         call r8_grafx4(theta,dist,distp,distsm,d_dist,ntk,
cdbg     >      'radians','m or cm',
cdbg     >      mnlab(1),mnlab(2),mnlab(3))
C
         do itk=1,ntk
            afac(itk)=distav/dist(itk)
            bfac(itk)=ONE
            gfac(itk)=distsm(itk)/dist(itk)
            pfac(itk)=afac(itk)*bfac(itk)*gfac(itk)
         enddo
C
cdbg         mnlab(2)=' '
cdbg         mnlab(3)='surface growth factors a,b,g, product'
cdbg         call r8_grafx4(theta,afac,bfac,gfac,pfac,ntk,
cdbg     >      ' ','m or cm',
cdbg     >      mnlab(1),mnlab(2),mnlab(3))
C
         do itk=1,ntk
            distp(itk)=dist(itk)
            dist(itk)=dist(itk)*(ONE+delta_x*pfac(itk))
            Rpts(itk)=Rcen+Rcoef(itk)*dist(itk)
            Zpts(itk)=Zcen+Zcoef(itk)*dist(itk)
         enddo
C
C  fft and then recompute distance
C
         call r8fftsc(Rpts, ntk, zs1d, zc1d, ierr)
         rmcx(jnext,0,1)=zc1d(1)/(2*ntk)
         do im=1,nmom
            rmcx(jnext,im,1) = zc1d(im+1)/ntk
            if(symflag) then
               rmcx(jnext,im,2) = ZERO
            else
               rmcx(jnext,im,2) = zs1d(im+1)/ntk
            endif
         enddo
 
         call r8fftsc(Zpts, ntk, zs1d, zc1d, ierr)
         if(symflag) then
            ymcx(jnext,0,1)=ZERO
         else
            ymcx(jnext,0,1)=zc1d(1)/(2*ntk)
         endif
         do im=1,nmom
            if(symflag) then
               ymcx(jnext,im,1) = ZERO
            else
               ymcx(jnext,im,1) = zc1d(im+1)/ntk
            endif
            ymcx(jnext,im,2) = zs1d(im+1)/ntk
         enddo
 
         call get_RZdist(jnext)         ! recompute from FFT expansion
 
      enddo
C
#ifdef __DEBUG
      mnlab=' '
      WRITE(MNLAB(1),8101)
 8101 FORMAT('r8bloata0 BLOAT SURFACE EXPANSION, final.')
      WRITE(MNLAB(2),8102) xiblo(jmax)
 8102 FORMAT('EXTRAPOLATION OF MHD SURFACE DATA, XImax=',1pe11.4)
      WRITE(MNLAB(3),8103)
 8103 FORMAT(1X,'(R,Y) BOUNDARY CONTOURS (MOMENTS EQUILIBRIUM)',14X)
      call r8_grfasym(rmcx,ymcx,mtbl,mimom,lcentr,jmax,nmom,
     >      'm or cm','m or cm',mnlab(1),mnlab(2),mnlab(3))
#endif
C
      RETURN
 
      contains
         !------------------------------------------------------------
         subroutine get_RZcen(jsurf)
 
            integer jsurf  ! surface where to get contour & distance profile
 
            Rcen=0; Zcen=0
            do itk=1,ntk
               call r8tkbmry(mtbl,jsurf,Theta(itk),rmcx,ymcx,mimom,nmom,
     >            Rpts(itk),Zpts(itk))
               Rcen=Rcen+Rpts(itk)
               Zcen=Zcen+Zpts(itk)
            enddo
            Rcen=Rcen/ntk
            Zcen=Zcen/ntk
 
         end subroutine get_RZcen
 
         !------------------------------------------------------------
         subroutine get_RZdist(jsurf)
 
            integer jsurf  ! surface where to get contour & distance profile
 
            do itk=1,ntk
               call r8tkbmry(mtbl,jsurf,Theta(itk),rmcx,ymcx,mimom,nmom,
     >            Rpts(itk),Zpts(itk))
               dist(itk)=sqrt((Rpts(itk)-Rcen)**2+(Zpts(itk)-Zcen)**2)
            enddo
 
         end subroutine get_RZdist
 
         !------------------------------------------------------------
         subroutine get_avdist
 
            real*8 zdth,zdthsum
            integer i
 
            distav=0
            d_distav=0
            zdthsum=0
 
            do i=1,ntk
               zdth=(ThetaA(i+1)-ThetaA(i-1))/2
               zdthsum=zdthsum+zdth
               distav=distav+dist(i)*zdth
               d_distav=d_distav+d_dist(i)*zdth
            enddo
            distav = distav/zdthsum
            d_distav = d_distav/zdthsum
 
         end subroutine get_avdist
 
         !------------------------------------------------------------
         subroutine find_distm1(jsurf)
 
            integer jsurf               ! surface where to get distances
 
c  find the intersection of the surface contour (jsurf) with the sequence
c  of line segments connecting the mag. axis to points on the boundary.
 
            real*8 rspl(4,ntk+1)
            real*8 zspl(4,ntk+1)
            real*8 wk(ntk+1)
 
            integer :: ibc = -1         ! periodic splines R(theta),Z(theta)
            real*8 :: zbc = 0.0_R8
 
            integer nn,ilinx
 
            integer, dimension(3) :: iselect = ( /1, 1, 0/ )
            real*8 rrv(3),zzv(3)        ! (R,dR/dtheta,-); (Z,dZ/dtheta,-)
 
            real*8 th                   ! (R(th),Z(th)) on line segment
            real*8 thp                  ! previous guess
            real*8 thn                  ! new guess, range not yet normalized.
 
            real*8 xR,xZ                ! vec perp {(Rbdy-Rcen),(Zbdy-Zcen)}
 
            real*8 f,dfdth
            !------------------------------------
 
            do itk=1,ntk
               call r8tkbmry(mtbl,jsurf,Theta(itk),rmcx,ymcx,mimom,nmom,
     >            Rspl(1,itk),Zspl(1,itk))
            enddo
 
            nn = ntk+1
            Rspl(1,nn)=Rspl(1,1)
            Zspl(1,nn)=Zspl(1,1)
 
            !  periodic splines R(theta), Z(theta)
            !  these describe the surface 1 zone width in from the
            !  boundary surface...
 
            call cspline(theta,nn,Rspl,ibc,zbc,ibc,zbc,wk,nn,ilinx,
     >         ierr)
            if(ierr.ne.0) call errmsg_exit('?cspline(Rspl) error.')
 
            call cspline(theta,nn,Zspl,ibc,zbc,ibc,zbc,wk,nn,ilinx,
     >         ierr)
            if(ierr.ne.0) call errmsg_exit('?cspline(Zspl) error.')
 
            do itk=1,ntk
               xR=-(Zpts(itk)-Zcen)
               xZ= (Rpts(itk)-Rcen)
 
               th = theta(itk)  ! initial guess (in range [theta(1),theta(nn)])
               thn = th         ! next guess (may be out of range)
               thp = c2pi*100   ! prev guess
               do
 
                  !  Newton search
 
                  call cspeval(th,iselect,rrv,theta,nn,ilinx,Rspl,
     >               ierr)
                  if(ierr.ne.0)
     >               call errmsg_exit('?cspeval(Rspl) error')
 
                  call cspeval(th,iselect,zzv,theta,nn,ilinx,Zspl,
     >               ierr)
                  if(ierr.ne.0)
     >               call errmsg_exit('?cspeval(Zspl) error')
 
                  if(abs(thn-thp).le.c2pi*1.0e-10_R8) exit
 
                  f=(rrv(1)-Rcen)*xR + (zzv(1)-Zcen)*xZ
                  dfdth = rrv(2)*xR + zzv(2)*xZ
 
                  thp=th
                  thn = th -f/dfdth
 
                  ! force new th in range, but save thn for comparision
 
                  if(thn.lt.theta(1)) then
                     th=thn+c2pi
                  else if(thn.gt.theta(nn)) then
                     th=thn-c2pi
                  else
                     th=thn
                  endif
 
               enddo
 
               distp(itk)=sqrt((Rcen-rrv(1))**2+(Zcen-zzv(1))**2)
            enddo
         end subroutine find_distm1
 
         !------------------------------------------------------------
         subroutine sm_dist
 
            real*8, dimension(:), allocatable :: dwk,tha
            integer iextnd,ilim1,ilim2,inum,i
 
            real*8, parameter :: zendc = 0.5_R8
 
            ! **quick** smooth of D(thetaA) data
 
            iextnd=ntk/2
            ilim1=1-iextnd
            ilim2=ntk+iextnd
            inum=ilim2-ilim1+1
 
            allocate(dwk(ilim1:ilim2),tha(ilim1:ilim2))
 
            tha(1:ntk)=thetaa(1:ntk)
            dwk(1:ntk)=dist(1:ntk)
            do i=1,iextnd
 
               tha(1-i)=tha(ntk+1-i)-c2pi
               tha(ntk+i)=tha(i)+c2pi
 
               dwk(1-i)=dwk(ntk+1-i)
               dwk(ntk+i)=dwk(i)
            enddo
 
            call r8_qksmooth(inum,tha,dwk, delsm, zendc,zendc, ierr)
            if(ierr.ne.0) call errmsg_exit('r8_qksmooth call error.')
 
            distsm(1:ntk)=dwk(1:ntk)
 
            deallocate(dwk,tha)
 
         end subroutine sm_dist
 
      END
