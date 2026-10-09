!******************** START FILE LOKJAC.FOR ; GROUP TKBLOAT ******************

!-- end FILE #BLOAT# ************************************
!FORMFEEDC-- START FILE #LOKJAC# ************************************
       LOGICAL FUNCTION LOKJAC(R0,DR0DXI,RJ,DRJDXI,YJ,DYJDXI,NMOM)
!
!
!       BOB MCCANN, 19-APR-85
!
!
!       --------
!       COMMENTS
!       --------
!
!
!	THIS CHECKS TO SEE if THE JACOBIAN D(R<Y)/D(XI,XI) CHANGES SIGN
!    AROUND THE GIVEN FLUX SURFACE.  LOKJAC=.T. if THE SIGN IS UNCHANGED,
!    AND LOKJAC=.F. if THE SIGN CHANGES
!    (THE TRANSFORMATION WENT SINGULAR).
!
!
!       -------------------
!       GLOBAL DECLARATIONS
!       -------------------
!
!
!***
!
!
!       ---------------
!       INPUT VARIABLES
!       ---------------
!
!
      integer ::	NMOM		!# OF MOMENTS OF R, Y
!
!    THE FOLLOWING REFER TO A GIVEN FLUX SURFACE:
!
      REAL	R0		!ZEROTH MOMENT
      REAL	DR0DXI		!D(R0)/DXI
      REAL	RJ(NMOM)	!J=1,NMOM MOMENTS OF R
      REAL	DRJDXI(NMOM)	!D(RJ)/DXI
!
      REAL	YJ(NMOM)	!K=1, NMOM MOMENTS OF Y
      REAL	DYJDXI(NMOM)	!D(YJ)/DXI
!
!
!       ----------------
!       OUTPUT VARIABLES
!       ----------------
!
!
!	LOGICAL LOKJAC	!=.T. if THE TRANSFORMATION IS OK
!			!=.F. if THE COORDINATES ARE SINGULAR
!
!
!       ------------------
!       LOCAL DECLARATIONS
!       ------------------
!
!
      parameter	(NTH = 127)	!# ANGLES TO USE IN FINDING EXTREMA
!
      REAL	ZC(NTH)		!COS(J*TH) ASSUMING NMOM.LE.NTH
      REAL	ZS(NTH)
!
      REAL	ZR77		!R(TH)
      REAL	ZY77		!Y(TH)
      REAL	ZDRDTH		!DR/DTH
      REAL	ZDYDTH		!DY/DTH
      REAL	ZDRDXI		!DR/DXI
      REAL	ZDYDXI		!DY/DXI
      REAL	ZJ		!JACOBIAN
!
!
!       ---------------
!       DATA STATEMENTS
!       ---------------
!
!
!***
!
!
!       -------------------
!       STATEMENT FUNCTIONS
!       -------------------
!
!
!***
!
!
!=========================================
!
!
!       ----------------------
!       0.1     INITIALIZATION
!       ----------------------
!
!
      ZPI = 4.0*ATAN(1.0)
!
!    LOOP OVER SOME ANGLES FROM 0 TO PI:
!
      do 69000 I = 1, NTH
!
        ZTH = ZPI*FLOAT(I-1)/FLOAT(NTH-1)
!
        do 100 J = 1, NMOM
          ZC(J) = COS(J*ZTH)
          ZS(J) = SIN(J*ZTH)
100     continue
!
!
!	----------------------------------------------
!	1.0	COMPUTE VARIOUS DERIVATIVES OF R AND Y
!	----------------------------------------------
!
!
10000 continue
!
!    HANDLE ZEROTH MOMENT AND INITIALIZE THE MOMENTS LOOP
!
      ZR77 = R0
      ZDRDXI = DR0DXI
      ZDRDTH = 0.0
!
      ZY77 = 0.0
      ZDYDTH = 0.0
      ZDYDXI = 0.0
!
!    ADD UP THE CONTRIBUTION FROM EACH MOMENT
!
      do 10100 J = 1, NMOM
!
        ZR77 = ZR77 + RJ(J)*ZC(J)
!
        ZDRDXI = ZDRDXI + DRJDXI(J)*ZC(J)
        ZDRDTH = ZDRDTH - RJ(J)*FLOAT(J)*ZS(J)
!
        ZY77 = ZY77 + YJ(J)*ZS(J)
!
        ZDYDXI = ZDYDXI + DYJDXI(J)*ZS(J)
        ZDYDTH = ZDYDTH + YJ(J)*FLOAT(J)*ZC(J)
!
10100 continue
!
!
!       -------------------------------------------------
!       2.0     COMPUTE THE JACOBIAN OF THIS FLUX SURFACE
!       -------------------------------------------------
!
!
20000   continue
!
      ZJ = ZDRDXI*ZDYDTH - ZDRDTH*ZDYDXI	!JACOBIAN
!
!
!	--------------------------------------------------
!	3.0	FIND THE MAXIMA AND MINIMA OF THE JACOBIAN
!	--------------------------------------------------
!
!
30000 continue
!
      if(I.EQ.1)	THEN	!INITIALIZE MAX AND MIN
!
        ZJMIN = ZJ
        THJMIN = ZTH
        ZJMAX = ZJ
        THJMAX = ZTH
!
      else			!FIND THE TRUE MIN AND MAX
!
        if(ZJ.LT.ZJMIN)	THEN
          ZJMIN = ZJ
          THJMIN = ZTH
        elseif(ZJ.GT.ZJMAX)	THEN
          ZJMAX = ZJ
          THJMAX = ZTH
        end if
!
      end if
!
69000 continue		!end OF ANGLE LOOP
!
!
!       ------------------------
!       7.0     CLEANUP & return
!       ------------------------
!
!
70000   continue
!
!
      if(ZJMIN*ZJMAX.GT.0.0)	THEN	!THEY HAVE THE SAME SIGN
          LOKJAC = .TRUE.
      else
          LOKJAC = .FALSE.
      end if
!
!
79000   continue        !ERROR returnS HERE
!
        return
!
!
!       ----------------------
!       8.0     ERROR HANDLING
!       ----------------------
!
!
!***
!
!
        end
!******************** end FILE LOKJAC.FOR ; GROUP TKBLOAT ******************
