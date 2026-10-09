C******************** START FILE LOKJAC.FOR ; GROUP TKBLOAT ******************
 
C-- END FILE #BLOAT# ************************************
CFORMFEEDC-- START FILE #LOKJAC# ************************************
       LOGICAL FUNCTION r8lokjac(R0,DR0DXI,RJ,DRJDXI,YJ,DYJDXI,NMOM)
C
C
C       BOB MCCANN, 19-APR-85
C
C
C       --------
C       COMMENTS
C       --------
C
C
C	THIS CHECKS TO SEE IF THE JACOBIAN D(R<Y)/D(XI,XI) CHANGES SIGN
C    AROUND THE GIVEN FLUX SURFACE.  LOKJAC=.T. IF THE SIGN IS UNCHANGED,
C    AND LOKJAC=.F. IF THE SIGN CHANGES
C    (THE TRANSFORMATION WENT SINGULAR).
C
C
C       -------------------
C       GLOBAL DECLARATIONS
C       -------------------
C
C
C***
C
C
C       ---------------
C       INPUT VARIABLES
C       ---------------
C
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      use iso_c_binding, only: fp => c_double 
      IMPLICIT NONE
      INTEGER nth,i,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 zpi,zth,zjmin,thjmin,zjmax,thjmax
!============
      INTEGER	NMOM		!# OF MOMENTS OF R, Y
C
C    THE FOLLOWING REFER TO A GIVEN FLUX SURFACE:
C
      REAL*8	R0		!ZEROTH MOMENT
      REAL*8	DR0DXI		!D(R0)/DXI
      REAL*8	RJ(NMOM)	!J=1,NMOM MOMENTS OF R
      REAL*8	DRJDXI(NMOM)	!D(RJ)/DXI
C
      REAL*8	YJ(NMOM)	!K=1, NMOM MOMENTS OF Y
      REAL*8	DYJDXI(NMOM)	!D(YJ)/DXI
C
C
C       ----------------
C       OUTPUT VARIABLES
C       ----------------
C
C
C	LOGICAL LOKJAC	!=.T. IF THE TRANSFORMATION IS OK
C			!=.F. IF THE COORDINATES ARE SINGULAR
C
C
C       ------------------
C       LOCAL DECLARATIONS
C       ------------------
C
C
      PARAMETER	(NTH = 127)	!# ANGLES TO USE IN FINDING EXTREMA
C
      REAL*8	ZC(NTH)		!COS(J*TH) ASSUMING NMOM.LE.NTH
      REAL*8	ZS(NTH)
C
      REAL*8	ZR77		!R(TH)
      REAL*8	ZY77		!Y(TH)
      REAL*8	ZDRDTH		!DR/DTH
      REAL*8	ZDYDTH		!DY/DTH
      REAL*8	ZDRDXI		!DR/DXI
      REAL*8	ZDYDXI		!DY/DXI
      REAL*8	ZJ		!JACOBIAN
C
C
C       ---------------
C       DATA STATEMENTS
C       ---------------
C
C
C***
C
C
C       -------------------
C       STATEMENT FUNCTIONS
C       -------------------
C
C
C***
C
C
C=========================================
C
C
C       ----------------------
C       0.1     INITIALIZATION
C       ----------------------
C
C
      ZPI = 4.0D0*ATAN(1.0D0)
C
C    LOOP OVER SOME ANGLES FROM 0 TO PI:
C
      DO 69000 I = 1, NTH
C
        ZTH = ZPI*REAL(I-1,kind=fp)/REAL(NTH-1,kind=fp)
C
        DO 100 J = 1, NMOM
          ZC(J) = COS(J*ZTH)
          ZS(J) = SIN(J*ZTH)
100     CONTINUE
C
C
C	----------------------------------------------
C	1.0	COMPUTE VARIOUS DERIVATIVES OF R AND Y
C	----------------------------------------------
C
C
10000 CONTINUE
C
C    HANDLE ZEROTH MOMENT AND INITIALIZE THE MOMENTS LOOP
C
      ZR77 = R0
      ZDRDXI = DR0DXI
      ZDRDTH = 0.0D0
C
      ZY77 = 0.0D0
      ZDYDTH = 0.0D0
      ZDYDXI = 0.0D0
C
C    ADD UP THE CONTRIBUTION FROM EACH MOMENT
C
      DO 10100 J = 1, NMOM
C
        ZR77 = ZR77 + RJ(J)*ZC(J)
C
        ZDRDXI = ZDRDXI + DRJDXI(J)*ZC(J)
        ZDRDTH = ZDRDTH - RJ(J)*REAL(J,kind=fp)*ZS(J)
C
        ZY77 = ZY77 + YJ(J)*ZS(J)
C
        ZDYDXI = ZDYDXI + DYJDXI(J)*ZS(J)
        ZDYDTH = ZDYDTH + YJ(J)*REAL(J,kind=fp)*ZC(J)
C
10100 CONTINUE
C
C
C       -------------------------------------------------
C       2.0     COMPUTE THE JACOBIAN OF THIS FLUX SURFACE
C       -------------------------------------------------
C
C
20000   CONTINUE
C
      ZJ = ZDRDXI*ZDYDTH - ZDRDTH*ZDYDXI	!JACOBIAN
C
C
C	--------------------------------------------------
C	3.0	FIND THE MAXIMA AND MINIMA OF THE JACOBIAN
C	--------------------------------------------------
C
C
30000 CONTINUE
C
      IF(I.EQ.1)	THEN	!INITIALIZE MAX AND MIN
C
        ZJMIN = ZJ
        THJMIN = ZTH
        ZJMAX = ZJ
        THJMAX = ZTH
C
      ELSE			!FIND THE TRUE MIN AND MAX
C
        IF(ZJ.LT.ZJMIN)	THEN
          ZJMIN = ZJ
          THJMIN = ZTH
        ELSEIF(ZJ.GT.ZJMAX)	THEN
          ZJMAX = ZJ
          THJMAX = ZTH
        ENDIF
C
      ENDIF
C
69000 CONTINUE		!END OF ANGLE LOOP
C
C
C       ------------------------
C       7.0     CLEANUP & RETURN
C       ------------------------
C
C
70000   CONTINUE
C
C
      IF(ZJMIN*ZJMAX.GT.0.0D0)	THEN	!THEY HAVE THE SAME SIGN
          r8lokjac = .TRUE.
      ELSE
          r8lokjac = .FALSE.
      ENDIF
C
C
79000   CONTINUE        !ERROR RETURNS HERE
C
        RETURN
C
C
C       ----------------------
C       8.0     ERROR HANDLING
C       ----------------------
C
C
C***
C
C
        END
C******************** END FILE LOKJAC.FOR ; GROUP TKBLOAT ******************
! 11Jan2003 fgtok -s r8_precision.sub r8bloat_names.sub "r8con.csh conversion"
! 11Jan2003 fgtok
