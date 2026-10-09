      subroutine r8bloats0(lcentr,nzones,xiblo,inbblo,mtbl,mj,
     >                  rmom0,rmom,ymom,nmom,mimom,nonlin,
     >                  rmb0,rmb2,ymb2)
C
C  bloat subroutine -- updown symmetric -- calling arguments, no COMMON!
C  dmc 27 Dec 1996
C
C  caller bloatx.for has some debug plotting stuff that may want to be
C  moved down into this routine in case the iteration loop fails...
C
C  this routine copies the inside-the-plasma moments set to the bloated
C  moments set, and uses a numerical algorithm to create well behaved
C  surfaces beyond the plasma bdy.  these surfaces provide a coordinate
C  system with reasonable properties (|J|>0 always) and a Fourier moments
C  representation, but the ones beyond the boundary are not real magnetic
C  surfaces!
C
C  when done ... set up spline coefficients for Fourier Spline -- used
C  by MOMRY
C
C input:
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
C  rmom0(..) -- 0th R moment set (inside plasma)
C  rmom(...) -- higher R moments set (inside plasma)
C  ymom(...) -- higher Y moments set (inside plasma)
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
C  rmb0(...) -- 0th R moments (extrapolated grid)
C  rmb2(...) -- higher R moments (extrapolated grid)
C  ymb2(...) -- higher Y moments (extrapolated grid)
C
C------------------------------------------------------------------
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      use iso_c_binding, only: fp => c_double
      IMPLICIT NONE
      INTEGER nzones,inbblo,mtbl,mj,nmom,mimom,nonlin,lcentr,jmax
      INTEGER inzp1,inzp2,lcm1,j,jzb,im,jbx0,iter,ibad,itermax,ierr
      INTEGER jb1,jb2,imxtra
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 ziter,zdx,zr0,zdr0dx
!============
      REAL*8 xiblo(mtbl)
C
      REAL*8 rmom0(mj)
      REAL*8 rmom(mimom,mj),ymom(mimom,mj)
C
      REAL*8 rmb0(mtbl)
      REAL*8 rmb2(mtbl,mimom)
      REAL*8 ymb2(mtbl,mimom)
C
C	--------------
      PARAMETER	(JMAX = 128)
C
C  UP-DOWN SYMMETRIC CASE :
C
      REAL*8	ZRJ(JMAX)               !R(J)
      REAL*8	ZYJ(JMAX)               !Y(J)
      REAL*8	ZDRDXI(JMAX)            !DR(J)/DXI
      REAL*8	ZDYDXI(JMAX)            !DY(J)/DXI
C
      LOGICAL	R8LOKJAC                  !LOGICAL FUNCTION - IS THE JACOBIAN OK?
C
C------------------------------------------------------------------
C
C   first copy the interior moments
C
C   NO. OF BDY'S INSIDE PLASMA; no. of zones PLUS 1
C
      INZP1=NZONES+1
      INZP2=NZONES+2
C
      LCM1=LCENTR-1
C
      DO 60 J = 1, INZP1
         JZB=J+LCM1
         RMB0(JZB) = RMOM0(JZB)
         DO 50 IM = 1, NMOM
            RMB2(JZB,IM) = RMOM(IM,JZB)
            YMB2(JZB,IM) = YMOM(IM,JZB)
 50      CONTINUE
 60   CONTINUE
C
C	DO THE EXTRAPOLATION INTO THE BLOATED ZONES
C
      JBX0=INZP2+LCM1
C
      iter=0
      ibad=0
      itermax=20
C
 99   continue
      ierr=0
      iter=iter+1
      ziter=REAL(iter,kind=fp)/REAL(itermax,kind=fp)
      if(ziter.gt.1.0D0) then
C
         WRITE(NONLIN,10150) XIBLO(JZB)
10150    FORMAT(' BLOAT ERROR: BAD BLOATED SURFACE AT XI=',1PE12.4
     1		,/,' THE JACOBIAN CHANGED SIGNS ON THE SURFACE'/
     1             ' Even on the final iteration (99 loop of BLOATX)')
         call errmsg_exit(' ... error in BLOATS0!')
C
      endif
C
      DO 10100 J = INZP2,INBBLO         !INTO GHOST ZONE
C
         JZB  = J + LCM1
         ZDX = XIBLO(JZB) - XIBLO(JZB-1)
C
C  DMC - CONSTRUCT AN EXTRAPOLATED SURFACE USING A LOCAL INTEGRATION
C   SCHEME.  MAKE THETA LINES SEEK A SIMPLE CIRCULAR SHAPE BUT ONLY
C   DIFFERENTIABLY, BASED ON DATA FROM THE PREVIOUS 2 SURFACES.
C
         CALL r8tksbcon(XIBLO,RMB0,MTBL,RMB2,YMB2,MIMOM,NMOM,
     >      JZB,JBX0,ZITER)
C
C  R8LOKJAC DATA
C
         JB1=JZB
         JB2=JZB-1
C
         ZR0=0.5D0*(RMB0(JB1)+RMB0(JB2))
C
         ZDR0DX=(RMB0(JB1)-RMB0(JB2))/ZDX
C
         DO 10090 IM=1,NMOM
C
            ZRJ(IM)=0.5D0*(RMB2(JB1,IM)+RMB2(JB2,IM))
            ZDRDXI(IM)=(RMB2(JB1,IM)-RMB2(JB2,IM))/ZDX
C
            ZYJ(IM)=0.5D0*(YMB2(JB1,IM)+YMB2(JB2,IM))
            ZDYDXI(IM)=(YMB2(JB1,IM)-YMB2(JB2,IM))/ZDX
C
10090    CONTINUE
C
         IMXTRA=NMOM
         IF(R8LOKJAC(ZR0,ZDR0DX,
     1      ZRJ,ZDRDXI,ZYJ,ZDYDXI,IMXTRA))	THEN
C
            GO TO 10100                 !COMPUTE THE NEXT SURFACE
C
         ELSE                           !THIS SURFACE WAS BAD, WARN & END LOOP
C
            ierr=1
            go to 9999
C
         ENDIF
C
10100 CONTINUE
C
 9999 CONTINUE
C
      if(ierr.eq.1) then
         ibad=iter
         go to 99                       ! try again to extrapolate surfaces
      endif
      if(iter.lt.itermax) then
         if(iter.eq.ibad+1) go to 99    ! try for some margin
      endif
C
      RETURN
      END
! 10Jan2003 fgtok -s r8_precision.sub r8bloat_names.sub "r8con.csh conversion"
! 10Jan2003 fgtok
