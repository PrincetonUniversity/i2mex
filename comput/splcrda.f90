!******************** START FILE SPLCRD.FOR ; GROUP TKSERVE ************
!.............................................................

        subroutine SPLCRDA (ZXI99, IJACI, &
                           LCENTR, MIMOM, NMOM, MTBL, XIBLO, &
                           RMCX,RMCXS,YMCX,YMCXS, &
      			   ZRMCX1,ZRMCX2,ZYMCX1,ZYMCX2, &
      		           ZRMCXP1,ZRMCXP2,ZYMCXP1,ZYMCXP2)
use iso_c_binding, only: fp => c_double

!
!  B.BALET - JET  MAR 1994 : MODifIED VERSION OF SPLCRD.FOR TO SUPPORT
!                            UP-doWN ASYMMETRIC CASE; SHOULD BE CALLED
!                            ONLY if LEVGEO = 6 or 7 or 8 or 9 (NLSYM = F)
!
!  dmc 24 Dec 1996:  this version does not depend on TRCOM:  all args passed.
!  ** updown asymmetric SPLCRDA **
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
!   RMCX(..) -- R moments array
!   RMCXS(..) -- R moments spline coeffs
!   YMCX(..) -- Y moments array
!   YMCXS(..) -- Y moments spline coeffs
!
!  output args:
!   ZRMCX1(..) -- R cos moments interpolated
!   ZRMCX2(..) -- R sin moments interpolated
!   ZYMCX1(..) -- Y cos moments interpolated
!   ZYMCX2(..) -- Y sin moments interpolated
!  (if IJAC set):  derivatives w.r.t. xi:
!   ZRMCXP1(..) -- R cos moments derivative
!   ZRMCXP2(..) -- R sin moments derivative
!   ZYMCXP1(..) -- Y cos moments derivative
!   ZYMCXP2(..) -- Y sin moments derivative
!
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

        REAL XIBLO(MTBL)
        REAL RMCX(MTBL,0:MIMOM,2),RMCXS(MTBL,0:MIMOM,2,3)
        REAL YMCX(MTBL,0:MIMOM,2),YMCXS(MTBL,0:MIMOM,2,3)

!   output arrays:

        DIMENSION ZRMCX1(0:MIMOM), ZRMCX2(0:MIMOM)
        DIMENSION ZYMCX1(0:MIMOM), ZYMCX2(0:MIMOM)
        DIMENSION ZRMCXP1(0:MIMOM), ZRMCXP2(0:MIMOM)
        DIMENSION ZYMCXP1(0:MIMOM), ZYMCXP2(0:MIMOM)

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
        if(IJAC.LE.1) THEN

!	EVALUATE SPLINE
!
           do 100 J = 0, NMOM

!	      RM(X) :
!
            ZRMCX1(J) = RMCX(ISL1,J,1) + ZDX99*(RMCXS(IS,J,1,1) + &
      		       ZDX99*(RMCXS(IS,J,1,2) + ZDX99*RMCXS(IS,J,1,3)))
              ZRMCX2(J) = RMCX(ISL1,J,2) + ZDX99*(RMCXS(IS,J,2,1) + &
                       ZDX99*(RMCXS(IS,J,2,2) + ZDX99*RMCXS(IS,J,2,3)))
!
!	      YM(X) :
!
              ZYMCX1(J) = YMCX(ISL1,J,1) + ZDX99*(YMCXS(IS,J,1,1) + &
                       ZDX99*(YMCXS(IS,J,1,2) + ZDX99*YMCXS(IS,J,1,3)))
              ZYMCX2(J) = YMCX(ISL1,J,2) + ZDX99*(YMCXS(IS,J,2,1) + &
                       ZDX99*(YMCXS(IS,J,2,2) + ZDX99*YMCXS(IS,J,2,3)))

 100       continue

        else
!
!  EVALUATE BY LINEAR INTERPOLATION
!
           ISP1=ISL1+1
           ZDX99I=ZDX99*ZDXII
!
           do 110 J = 0, NMOM

              ZRMCX1(J)= RMCX(ISL1,J,1) + &
                 ZDX99I*(RMCX(ISP1,J,1)-RMCX(ISL1,J,1))
              ZRMCX2(J)= RMCX(ISL1,J,2) + &
                 ZDX99I*(RMCX(ISP1,J,2)-RMCX(ISL1,J,2))
              ZYMCX1(J)= YMCX(ISL1,J,1) + &
                 ZDX99I*(YMCX(ISP1,J,1)-YMCX(ISL1,J,1))
              ZYMCX2(J)= YMCX(ISL1,J,2) + &
                 ZDX99I*(YMCX(ISP1,J,2)-YMCX(ISL1,J,2))

 110       continue

        end if

!..................................

        if (IJAC.EQ.0) return
        if (IJAC.EQ.2) return

!..................................

        if(IJAC.LE.1) THEN
!
!  SPLINES
!
           do 200 J = 0, NMOM

!     D[RM(X)]/DX
              ZRMCXP1(J) = RMCXS(IS,J,1,1) + ZDX99*(2.*RMCXS(IS,J,1,2) + &
                 3.*ZDX99*RMCXS(IS,J,1,3))
              ZRMCXP2(J) = RMCXS(IS,J,2,1) + ZDX99*(2.*RMCXS(IS,J,2,2) + &
                 3.*ZDX99*RMCXS(IS,J,2,3))

!     D[YM(X)]/DX
              ZYMCXP1(J) = YMCXS(IS,J,1,1) + ZDX99*(2.*YMCXS(IS,J,1,2) + &
                 3.*ZDX99*YMCXS(IS,J,1,3))
              ZYMCXP2(J) = YMCXS(IS,J,2,1) + ZDX99*(2.*YMCXS(IS,J,2,2) + &
                 3.*ZDX99*YMCXS(IS,J,2,3))

 200       continue

        else
!
!  LINEAR INTERPOLATION
!
           do 210 J = 0, NMOM

              ZRMCXP1(J) = (RMCX(ISP1,J,1)-RMCX(ISL1,J,1))*ZDXII
              ZRMCXP2(J) = (RMCX(ISP1,J,2)-RMCX(ISL1,J,2))*ZDXII
              ZYMCXP1(J) = (YMCX(ISP1,J,1)-YMCX(ISL1,J,1))*ZDXII
              ZYMCXP2(J) = (YMCX(ISP1,J,2)-YMCX(ISL1,J,2))*ZDXII

 210       continue

        end if
!
        return
        end
!******************** end FILE SPLCRD.FOR ; GROUP TKSERVE **************
