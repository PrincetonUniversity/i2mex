C-----------------------------------------------------------------
C  TDB_PRFNXP
C   FIND LEAST POINT IN PROFILE THAT IS GREATER THAN PASSED VALUE
C   ALSO REQUIRE POINT TO BE IN REGION OF CURVE WITH POSITIVE SLOPE,
C   OR THE CURVES ABSOLUTE MAXIMUM.
C
      SUBROUTINE TDB_PRFNXP(ZY,ZX,INUM,ZBASE,ZX1,ILOC,ZTOL)
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER inum,iloc,imax,isets,i,ilm1,ilp1,ifnd,ix,imin
      integer lunmsg_tdb
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 zx,zbase,zx1,ztol,zy,zcd,zc0,zc1,zx0,zy0,zmax,zmin
!============
      DIMENSION ZX(INUM),ZY(INUM)
C
C  ZX,ZY (INUM)  INPUT PROFILE
C  ZBASE  INPUT BASE VALUE
C  ILOC   OUTPUT INDEX LOCATION
C  ZX1    OUTPUT FIRST X LOCATION
C
C  FIND LEAST POINT ZY(ILOC) .GT. ZBASE
C  RETURN ILOC=0 IF NO SUCH POINT EXISTS
C
      ZCD=.5_R8*ABS(ZX(INUM)-ZX(1))
      ZC0=ZX(1)-2._R8*ZCD
      ZC1=ZC0+ZCD
 1    CONTINUE
C  LOCATE ALSO FIRST RADIUS WHERE PROFILE SLOPE IS INCREASING
      ZX0=ZC0
      zx1=zx0
      zy0=0.0_R8
      ILOC=0
      ZMAX=0.0_R8
      IMAX=0
      IMIN=0
      isets=0
      DO 10 I=1,INUM
         if(i.eq.1) then
            imin=1
            zmin=zy(1)
         else if(zy(i).lt.zmin) then
            imin=i
            zmin=zy(i)
         endif
C  INCR. SLOPE?
         IF(isets.gt.0) GO TO 2
C
         IF(I.EQ.INUM) THEN
            ZX0=ZC1
            zy0=zmin
         else
            IF(ZY(I+1).LE.ZY(I)) GO TO 2
            ZX0=ZC1
            IF(I.GT.1) ZX0=.5_R8*(ZX(I-1)+ZX(I))
            ZY0=ZY(I)
            isets=1
         endif
C
 2       CONTINUE
         IF(ZY(I).LE.ZMAX) GO TO 4
         ZMAX=ZY(I)
         IMAX=I
C
 4       CONTINUE
         IF(ZY(I).LE.ZBASE) GO TO 10
         IF(ILOC.GT.0) GO TO 6
         ILOC=I
C
 6       CONTINUE
         IF(ZY(I).LT.ZY(ILOC)) ILOC=I
 10   CONTINUE
C
      IF(ILOC.EQ.0) GOTO 1000
      IF(IPRFEQ(ZY(ILOC),ZY(IMAX),ZTOL)) ILOC=IMAX
      ZX1=ZX(ILOC)
      ZBASE=ZY(ILOC)*(1.0_R8+ZTOL)
C
      IF(ILOC.EQ.IMAX) GOTO 1000
C---
C
      ILM1=max(1,(ILOC-1))
      IF((ILOC.EQ.1).AND.(ZY(ILOC+1).GT.ZY(ILOC))) GOTO 1000
      IF(ILOC.EQ.INUM) GO TO 100
      IF(ZY(ILOC).EQ.ZY(ILOC+1)) ILOC=ILOC+1
      ILP1=min(INUM,(ILOC+1))
C  MAKE SURE POINT IS NOT A LOCAL SINGULARITY
      IF((ZY(ILOC).LT.ZY(ILM1)).AND.(ZY(ILOC).LT.ZY(ILP1)))
     >     GO TO 1
      IF((ZY(ILOC).GT.ZY(ILM1)).AND.(ZY(ILOC).GT.ZY(ILP1)))
     >     GO TO 1
C  BE SURE TO LOCATE CORRECT FIRST RADIAL LOCATION--
C  IN CASE OF HOLLOW PROFILE
 100  CONTINUE
      IF(ZY(ILOC).LE.ZY0) GO TO 1
      CALL TDB_PRFNXX(ZY,ZX,INUM,ZY(ILOC),ZX0,ZX1,IFND)
      GOTO 1000
C-----
 1000 CONTINUE
      RETURN
C-----

      contains
        logical function iprfeq(zt1,zt2,ztol)
        real*8, intent(in) :: zt1,zt2,ztol

        real*8 :: z1,z2

        Z1=min(ZT1,ZT2)
        Z2=max(ZT1,ZT2)
        IPRFEQ=.FALSE.
        IF(Z2.LE.(1.0_R8+ZTOL)*Z1) IPRFEQ=.TRUE.

        end function iprfeq

      END
