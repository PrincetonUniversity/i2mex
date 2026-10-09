C******************** START FILE PLCMJR.FOR ; GROUP RPLOT_SUB  ******************
C==============================================================
C  PLCMJR   Change X axis in RPLOT calculator to Major radius.
 
 
 
      SUBROUTINE PLCMJR (Jtype,ier)
 
 
 
C       Updates:
C       DMC 9/15/99  added IER argument; implement as %RMJMAP command.
C       TBT 4/14/94  Commented out IPX and IPXMR which aren't used.
C       tbt 3/08/94  Changed algorithm for center point when data
C                    is zone centered.
C       TBT 9/29/92  Made for TPLT3D.for
C
C	Argument: Input: Jtype - X axis type of accumulator in calculator.
C                                Must be 1 or 2.
C                                1 = Zone center
C                                2 = Zone Boundary.
C                                5 = Major radius.
C                 Output: Jtype = 5, or whatever RMAJM type is.
C
C  COMMON BLOCKS---
      use datmgr_mod
      use cplotr_mod
      use rpcalc_mod
C
C  LOCAL--
C
      CHARACTER*1 IANS
C
      DIMENSION ZWK(NR0)
 
      Real DatTemp(NR0)        ! Temporary buffer for X axis data.
      Real Data1, Data2        ! First and second points to calculate center
 
C
C----------------------------------------------------------------
C
C  START OF EXECUTABLE CODE--
C
      ier = 0
      If (Jtype .Eq. 0)  Then
          Call ZERMSG(' PLCMJR:"$" is UNDEFINED')
          ier=ier+1
      End if      ! Jtype
 
	
      If (Jtype .Ne. 1 .AND. Jtype .Ne. 2)  Then
          Call ZERMSG(' PLCMJR: "$" is not vs. flux zone ctr/bdy')
          ier=ier+1
      End if      ! Jtype
 
 
C   "RMAJM" MUST BE AMONG THE COLLECTION OF X AXES
      IF(.NOT.LRMAJM) THEN
         CALL ZERMSG(' MAJOR RADIUS PLOT DATA "RMAJM" NOT FOUND')
         ier=ier+1
      ENDIF
      INXMR=NZONEX(NTYPMR)
      IF((NZONEX(1).NE.NZONEX(2)).OR.(INXMR.NE.(2*NZONEX(1)+1)))THEN
         CALL ZERMSG(' MAJOR RADIUS "RMAJM" HAS WRONG NO. OR PTS.')
         ier=ier+1
      ENDIF
C
      if(ier.gt.0) return
C
      JTYP = Jtype
C
C
C  READ THE DATA AND X AXIS IF APPLICABLE
C  X AXIS
      IND1=0
      IND2=0
      IF(NLXFOT(JTYP)) THEN
        CALL DMGXOT(JTYP,IND1,IND2)
        CALL PLPRIO(6,IND1)
        CALL PLPRIO(6,IND2)
      ENDIF
C
C  RMAJM  --  OR ALTERNATE X AXIS SUBSTITUTED FOR RMAJM...
        INDMR=0
        CALL DMGXOT(NTYPMR,INDMR,IDUM)
        CALL PLPRIO(6,INDMR)
 
C  RESTORE STANDARD PRIORITY FOR FCNS READ
      CALL PLPRIO(5,IND1)
      CALL PLPRIO(5,IND2)
      CALL PLPRIO(5,INDMR)
C
C
        CALL DMDLOC('%WRK1',IND1,ISIZ1,IWRK1)
      CALL DMDLOC('%WRK2',IND2,ISIZ2,IWRK2)
C
      INX=NZONEX(JTYP)
      IF(.NOT.NLXFOT(JTYP)) CALL copyr4(XARRY(1,JTYP),XF,INX)
 
C
 
C	.Copy accumulator data from Workspace 2 to workspace 1.
            Call copyr4(DatBuf(Iwrk2), DatBuf(Iwrk1), INX*NTR)
 
      IPF = Iwrk1      ! tbt
 
      DO 180 IT=1,NTR
 
       IPFL=IPF+(IT-1)*INX
C        CHECK
C        ... for ZONE BDYS IF MAPPING TO RMAJM GRID
       IF(2 .EQ. JTYP) THEN
C         .Accumulator is zone boundary - just copy data.
        CALL copyr4(DATBUF(IPFL),DatTemp(1),INX)
 
       ELSE
C         .Switch accumulator from zone center to zone boundary.
        CALL XINTZB(DATBUF(IPFL),DatTemp(1),INX)
       ENDIF
 
 
C         MAP TO MAJ. RADIUS GRID
 
C         AXIS EXTRAPOLATION
 
C	 .Fit data to parabola (f(x) = a*x*x + b) for center point
       IF(2 .EQ. JTYP) THEN
C         .Accumulator is zone boundary - find center point.
        Data1 = DatTemp(1)
        Data2 = DatTemp(2)
        ZData0 = Data1 + 0.3333333*(Data1-Data2)  ! Full grid box
 
       ELSE
C         .Find center point for zone center.
        Data1 = DatBuf(IPFL)
        Data2 = DatBuf(IPFL+1)
        ZData0 = Data1 + 0.125*(Data1-Data2)      ! Half grid box.
       ENDIF
 
 
        IPWL=IWRK2+(IT-1)*INXMR
        ICEN=IPWL+INX
        DATBUF(ICEN)=ZDATA0
 
        DO 188 IX=1,INX
          DATBUF(ICEN-IX)=DatTemp(IX)
          DATBUF(ICEN+IX)=DatTemp(IX)
 188    CONTINUE
 
C        LOOP END
 180  CONTINUE
C
C  HAVING REMAPPED DATA, PASS RMAJM AS X AXIS TO Calculator
C
 
      Jtype = Ntypmr
C
      Return
C
      	END
C
C******************** END FILE PLCMJR.FOR ; GROUP RPLOT_SUB  ******************
 
