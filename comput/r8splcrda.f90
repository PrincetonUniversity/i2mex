subroutine r8splcrda (ZXI99, IJACI, &
     LCENTR, MIMOM, NMOM, MTBL, XIBLO, &
     RMCX,RMCXS,YMCX,YMCXS, &
     ZRMCX1,ZRMCX2,ZYMCX1,ZYMCX2, &
     ZRMCXP1,ZRMCXP2,ZYMCXP1,ZYMCXP2)
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
  !	COORD "ZXI99", RETURN THE MHD FOURIER COEFFICIENTS (AND THEIR
  !	DERIVATIVES) DESCRIBING THE PLASMA SURFACE ON THAT COORD.
  !
  !  MOD DMC APRIL 1990 -- LINEAR INTERPOLATION OPTION
  !
  !  if IJAC=0 OR IJAC=1 USE THE SPLINES
  !  if IJAC=2 OR IJAC=3 do LINEAR INTERPOLATION INSTEAD
  !
  !  if IJAC=1 OR IJAC=3 EVALUATE DATA FOR THE JACOBIAN
  !
  use iso_c_binding, only: fp => c_double
  implicit none
  !   input arrays:
  integer :: ijaci,lcentr,mimom,nmom,mtbl,ijac,is,isl1,j,isp1
  real(fp) :: zxi99,zdxii,zdx99,zdx99i
  real(fp) :: XIBLO(MTBL)
  real(fp) :: RMCX(MTBL,0:MIMOM,2),RMCXS(MTBL,0:MIMOM,2,3)
  real(fp) :: YMCX(MTBL,0:MIMOM,2),YMCXS(MTBL,0:MIMOM,2,3)

  !   output arrays:
  real(fp), dimension(0:mimom) :: zrmcx1,zrmcx2
  real(fp), dimension(0:mimom) :: zymcx1,zymcx2
  real(fp), dimension(0:mimom) :: zrmcxp1,zrmcxp2
  real(fp), dimension(0:mimom) :: zymcxp1,zymcxp2
  !
  !-------------------------------------
  !
  IJAC=abs(IJACI)
  !
  ZDXII= 1.0D0/ (XIBLO(LCENTR+1)-XIBLO(LCENTR))
  IS = ZXI99 * ZDXII + 1
  !
20 continue
  ISL1=IS+LCENTR-1
  if(zxi99.le.xiblo(isl1)) then
    is=is-1
    go to 20
  end if
  if(zxi99.gt.xiblo(isl1+1)) then
    is=is+1
    go to 20
  end if
  !
  ZDX99 = ZXI99 - XIBLO(ISL1)
  !
  if(IJAC.LE.1) THEN
 
    !	EVALUATE SPLINE
    !
    do J = 0, NMOM

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

    end do

  else
    !
    !  EVALUATE BY LINEAR INTERPOLATION
    !
    ISP1=ISL1+1
    ZDX99I=ZDX99*ZDXII
    !
    do J = 0, NMOM

      ZRMCX1(J)= RMCX(ISL1,J,1) + &
           ZDX99I*(RMCX(ISP1,J,1)-RMCX(ISL1,J,1))
      ZRMCX2(J)= RMCX(ISL1,J,2) + &
           ZDX99I*(RMCX(ISP1,J,2)-RMCX(ISL1,J,2))
      ZYMCX1(J)= YMCX(ISL1,J,1) + &
           ZDX99I*(YMCX(ISP1,J,1)-YMCX(ISL1,J,1))
      ZYMCX2(J)= YMCX(ISL1,J,2) + &
           ZDX99I*(YMCX(ISP1,J,2)-YMCX(ISL1,J,2))

    end do

  end if
 
  !..................................
  if (IJAC.EQ.0) RETURN
  if (IJAC.EQ.2) RETURN
  !..................................

  if(IJAC.LE.1) THEN
    !
    !  SPLINES
    !
    do J = 0, NMOM

      !     D[RM(X)]/DX
      ZRMCXP1(J) = RMCXS(IS,J,1,1) + &
           ZDX99*(2._fp*RMCXS(IS,J,1,2) + &
           3._fp*ZDX99*RMCXS(IS,J,1,3))
      ZRMCXP2(J) = RMCXS(IS,J,2,1) + &
           ZDX99*(2._fp*RMCXS(IS,J,2,2) + &
           3._fp*ZDX99*RMCXS(IS,J,2,3))
 
      !     D[YM(X)]/DX
      ZYMCXP1(J) = YMCXS(IS,J,1,1) + &
           ZDX99*(2._fp*YMCXS(IS,J,1,2) + &
           3._fp*ZDX99*YMCXS(IS,J,1,3))
      ZYMCXP2(J) = YMCXS(IS,J,2,1) + &
           ZDX99*(2._fp*YMCXS(IS,J,2,2) + &
           3._fp*ZDX99*YMCXS(IS,J,2,3))
 
    end do
 
  else
    !
    !  LINEAR INTERPOLATION
    !
    do J = 0, NMOM
 
      ZRMCXP1(J) = (RMCX(ISP1,J,1)-RMCX(ISL1,J,1))*ZDXII
      ZRMCXP2(J) = (RMCX(ISP1,J,2)-RMCX(ISL1,J,2))*ZDXII
      ZYMCXP1(J) = (YMCX(ISP1,J,1)-YMCX(ISL1,J,1))*ZDXII
      ZYMCXP2(J) = (YMCX(ISP1,J,2)-YMCX(ISL1,J,2))*ZDXII

    end do

  end if
  !
  RETURN
END subroutine r8splcrda
