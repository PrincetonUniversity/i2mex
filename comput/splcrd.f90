!******************** START FILE SPLCRD.FOR ; GROUP TKSERVE ************
!.............................................................

        subroutine SPLCRD (ZXI99, IJACI, &
                           LCENTR, MIMOM, NMOM, MTBL, XIBLO, &
                           RMB0,RMB0S,RMB2,RMB2S,YMB2,YMB2S, &
      			   ZRMB0X, ZRMB2X, ZYMB2X, &
      		           ZRMB0P, ZRMB2P, ZYMB2P)
use iso_c_binding, only: fp => c_double

!
!  dmc 24 Dec 1996:  this version does not depend on TRCOM
!  ** updown symmetric SPLCRD **
!
!  call from MOMRY/MOMRZ; caller does arg error checking!
!
!  input args:
!
!   ZXI99 -- "xi" (radial coordinate) **** must be in range ****
!     ==> interpolate moments set to "xi"
!   IJACI -- IJAC flag as input (see below)
!
!   LCENTR -- 1 + index offset in XIBLO, RMB0, RMB2, YMB2, rel. to.
!      RMB0S,RMB2S,YMB2S
!      (for TRANSP compatibility; set to 1 if all arrays are on
!      the same indexing grid)
!
!   MIMOM -- max no. of moments (array dimension)
!   NMOM  -- actual no. of moments in use
!   MTBL --- max no. of radial points including "bloat" (array dimension)
!
!   XIBLO(..) -- "xi" grid points, including bloat region beyond plasma bdy
!
!   RMB0(..) -- 0'th R moment array
!   RMB0S(..) -- 0'th R moment spline coeffs
!   RMB2(..) -- higher R moments array
!   RMB2S(..) -- higher R moment spline coeffs
!   YMB2(..) -- higher Y moments array
!   YMB2S(..) -- higher Y moment spline coeffs
!
!  output args:
!   ZRMB0X -- interpolated 0'th R moment
!   ZRMB2X(..) -- interpolated higher R moments
!   ZYMB2X(..) -- interpolated higher Y moments
!  (if IJAC set):
!   ZRMB0P -- 0'th R moment dR0/dxi
!   ZRMB2P(..) -- higher R moments derivatives w.r.t. xi
!   ZYMB2P(..) -- higher Y moments derivatives w.r.t. xi
!
!	EVALUATE SPLINE COEFFICIENTS DIRECTLY ... GIVEN A RADIAL
!	COORD "ZXI99", return THE MHD FOURIER COEFFICIENTS (AND THEIR
!	DERIVATIVES) DESCRIBING THE PLASMA SURFACE ON THAT COORD.
!
!  MOD DMC APRIL 1990 -- LINEAR INTERPOLATION OPTION
!
!  if IJAC=0 OR IJAC=1 USE THE SPLINES
!  if IJAC=2 OR IJAC=3 do LINEAR INTERPOLATION INSTEAD
!
!  if IJAC=1 OR IJAC=3 EVALUATE DATA FOR THE JACOBIAN
!
!   input arrays:

        REAL RMB0(MTBL),RMB0S(MTBL,3),XIBLO(MTBL)
        REAL RMB2(MTBL,MIMOM),RMB2S(MTBL,MIMOM,3)
        REAL YMB2(MTBL,MIMOM),YMB2S(MTBL,MIMOM,3)

!   output arrays:

        DIMENSION ZRMB2X(MIMOM), ZYMB2X(MIMOM)
        DIMENSION ZRMB2P(MIMOM), ZYMB2P(MIMOM)

!
!-------------------------------------
!
        IJAC=abs(IJACI)
!
        ZDXII= 1.0 / (XIBLO(LCENTR+1)-XIBLO(LCENTR))
        IS = ZXI99 * ZDXII + 1
!
 20     continue
        ISL1=IS+LCENTR-1
        if(zxi99.le.xiblo(isl1)) then
           is=is-1
           goto 20
        end if
        if(zxi99.gt.xiblo(isl1+1)) then
           is=is+1
           goto 20
        end if
!
        ZDX99 = ZXI99 - XIBLO(ISL1)
!
!  EVALUATE MOMENTS INTERPOLATION
!
        if(IJAC.LE.1) THEN

!	EVALUATE SPLINE
!
!	  R0(X)
!
           ZRMB0X = RMB0(ISL1) + ZDX99*(RMB0S(IS,1) + &
      			    ZDX99*(RMB0S(IS,2) + &
      			    ZDX99*RMB0S(IS,3)))


           do 100 J = 1, NMOM

!	      RM(X)
            ZRMB2X(J) = RMB2(ISL1,J) + ZDX99*(RMB2S(IS,J,1) + &
      			   ZDX99*(RMB2S(IS,J,2) + ZDX99*RMB2S(IS,J,3)))
!	      YM(X)
            ZYMB2X(J) = YMB2(ISL1,J) + ZDX99*(YMB2S(IS,J,1) + &
      			   ZDX99*(YMB2S(IS,J,2) + ZDX99*YMB2S(IS,J,3)))

 100       continue

        else
!
!  EVALUATE BY LINEAR INTERPOLATION
!
           ISP1=ISL1+1
           ZDX99I=ZDX99*ZDXII
!
           ZRMB0X=RMB0(ISL1)+ZDX99I*(RMB0(ISP1)-RMB0(ISL1))
!
           do 110 J = 1, NMOM

            ZRMB2X(J)= RMB2(ISL1,J)+ZDX99I*(RMB2(ISP1,J)-RMB2(ISL1,J))
            ZYMB2X(J)= YMB2(ISL1,J)+ZDX99I*(YMB2(ISP1,J)-YMB2(ISL1,J))

 110       continue

        end if

!..................................

      if (IJAC.EQ.0) return
      if (IJAC.EQ.2) return

!..................................

!
!  EVALUATE DERIVATIVES
!

        if(IJAC.LE.1) THEN
!
!  SPLINES
!
!	  D[R0(X)]/DX
           ZRMB0P = RMB0S(IS,1) + ZDX99*(2.*RMB0S(IS,2) + &
      		     3.*ZDX99*RMB0S(IS,3))

           do 200 J = 1, NMOM

!	  D[RM(X)]/DX
            ZRMB2P(J) = RMB2S(IS,J,1) + ZDX99*(2.*RMB2S(IS,J,2) + &
      			   3.*ZDX99*RMB2S(IS,J,3))
!	  D[YM(X)]/DX
            ZYMB2P(J) = YMB2S(IS,J,1) + ZDX99*(2.*YMB2S(IS,J,2) + &
      			   3.*ZDX99*YMB2S(IS,J,3))

 200       continue

        else
!
!  LINEAR INTERPOLATION
!
           ZRMB0P=(RMB0(ISP1)-RMB0(ISL1))*ZDXII
!
           do 210 J = 1, NMOM

            ZRMB2P(J) = (RMB2(ISP1,J)-RMB2(ISL1,J))*ZDXII
            ZYMB2P(J) = (YMB2(ISP1,J)-YMB2(ISL1,J))*ZDXII

 210       continue

        end if
!
        return
        end
!******************** end FILE SPLCRD.FOR ; GROUP TKSERVE **************
