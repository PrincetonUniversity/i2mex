       LOGICAL FUNCTION LOKJACA &
                        (RJC,RJS,DRCDXI,DRSDXI,YJC,YJS,DYCDXI,DYSDXI,NM)
!
!
!       B. BALET, JAN 1994 : MODifICATIONS TO THE FUNCTION LOKJAC WRITTEN BY
!                            BOB MCCANN FOR THE UP-doWN ASYMMETRIC CASE
!
!       --------
!       COMMENTS
!       --------
!
!
!	THIS CHECKS TO SEE if THE JACOBIAN D(R<Y)/D(XI,XI) CHANGES SIGN
!    AROUND THE GIVEN FLUX SURFACE.  LOKJACA=.T. if THE SIGN IS UNCHANGED,
!    AND LOKJACA=.F. if THE SIGN CHANGES
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
      integer ::	NM  		!# OF MOMENTS OF R, Y
!
!    THE FOLLOWING REFER TO A GIVEN FLUX SURFACE:
!
      REAL	RJC(0:NM) 	!J=0,NM MOMENTS OF R : COS TERMS
        REAL	RJS(0:NM) 	!J=0,NM MOMENTS OF R : SIN TERMS
      REAL	DRCDXI(0:NM)	!D(RJC)/DXI
      REAL	DRSDXI(0:NM)	!D(RJS)/DXI
!
        REAL    YJC(0:NM)         !J=0,NM MOMENTS OF Y : COS TERMS
        REAL    YJS(0:NM)         !J=0,NM MOMENTS OF Y : SIN TERMS
        REAL    DYCDXI(0:NM)      !D(YJC)/DXI
        REAL    DYSDXI(0:NM)      !D(YJS)/DXI
!
!
!       ----------------
!       OUTPUT VARIABLES
!       ----------------
!
!
!	LOGICAL LOKJACA	!=.T. if THE TRANSFORMATION IS OK
!			!=.F. if THE COORDINATES ARE SINGULAR
!
!
!       ------------------
!       LOCAL DECLARATIONS
!       ------------------
!
!
      parameter	(NTH = 254)	!# ANGLES TO USE IN FINDING EXTREMA
!
      REAL	ZC(0:NTH)		!COS(J*TH) ASSUMING NM.LE.NTH
      REAL	ZS(0:NTH)         !SIN(J*TH) ASSUMING NM.LE.NTH
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
      ZTWOPI = 8.0*ATAN(1.0)
!
!    LOOP OVER SOME ANGLES FROM 0 TO TWO*PI:
!
      do 69000 I = 1, NTH
!
        ZTH = ZTWOPI*FLOAT(I-1)/FLOAT(NTH-1)
!
        do 100 J = 0, NM
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
!    INITIALIZE THE MOMENTS LOOP
!
      ZR77 = 0.0
      ZDRDXI = 0.0
      ZDRDTH = 0.0
!
      ZY77 = 0.0
      ZDYDTH = 0.0
      ZDYDXI = 0.0
!
!    ADD UP THE CONTRIBUTION FROM EACH MOMENT
!
      do 10100 J = 0, NM
!
        ZR77 = ZR77 + RJC(J)*ZC(J) + RJS(J)*ZS(J)
!
        ZDRDXI = ZDRDXI + DRCDXI(J)*ZC(J) + DRSDXI(J)*ZS(J)
        ZDRDTH = &
          ZDRDTH - RJC(J)*FLOAT(J)*ZS(J) + RJS(J)*FLOAT(J)*ZC(J)
!
          ZY77 = ZY77 + YJC(J)*ZC(J) + YJS(J)*ZS(J)
!
          ZDYDXI = ZDYDXI + DYCDXI(J)*ZC(J) + DYSDXI(J)*ZS(J)
          ZDYDTH = &
          ZDYDTH - YJC(J)*FLOAT(J)*ZS(J) + YJS(J)*FLOAT(J)*ZC(J)
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
          LOKJACA = .TRUE.
      else
          LOKJACA = .FALSE.
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
