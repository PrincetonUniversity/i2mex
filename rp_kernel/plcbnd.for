C******************** START FILE PLCBND.FOR ; GROUP RPLOT_SUB  ******************
C==============================================================
C  PLCBND   Change X axis in RPLOT calculator to zone boundary
C                                    from Major radius.
 
      SUBROUTINE PLCBND (Jtype, itarg, imap, ier)
 
C       Updates:
C       dmc  9/15/99 added IMAP and IER arguments; removed UREAD call.
C       TBT 12/14/94 Created from PLCMJR.for
C
C
C	Argument: Input: Jtype - X axis type of accumulator in calculator.
C                                Must be NTYPMR.
C                                1 means x is Zone center
C                                2 means x is Zone Boundary.
C                                  ...
C                                NTYPMR means x is Major radius.
C                 Output: Jtype = itarg (if successful)
C
C                 Input:  Itarg=1:  output is to be zone ctr'd
C                         Itarg=2:  output is to be bdy ctr'd
C
C                 Input:  Imap=1:  map R>Raxis data to xb
C                         Imap=2:  map R<Raxis data to xb
C                         Imap=3:  map by averaging both sides
C
C                 Output:  ier, completion code, 0 = ok
C
C  COMMON BLOCKS---
      use datmgr_mod
      use cplotr_mod
      use rpcalc_mod
C
C  LOCAL--
C
C
C----------------------------------------------------------------
C
C  START OF EXECUTABLE CODE--
C
C   "RMAJM" MUST BE AMONG THE COLLECTION OF X AXES
      ier = 0
C
      IF(.NOT.LRmajm) THEN
         CALL ZERMSG(' MAJOR RADIUS PLOT DATA "RMAJM" NOT FOUND.')
         Ier=Ier+1
      ENDIF
      INXMR=NZONEX(NTYPMR)
      IF((NZONEX(1).NE.NZONEX(2)).OR.(INXMR.NE.(2*NZONEX(1)+1)))THEN
         CALL ZERMSG(' MAJOR RADIUS "RMAJM" HAS WRONG NO. OR PTS.')
         Ier=Ier+1
      ENDIF
	
      If (Jtype .Eq. 0)  Then
         Call ZERMSG(' PLCBND:"$" is UNDEFINED.')
         Ier=Ier+1
      End if                            ! Jtype
 
	
      If (Jtype .Ne. Ntypmr)  Then
         Call ZERMSG(' PLCBND: "$" is not vs. Major Radius.')
         Ier=Ier+1
      End if                            ! Jtype
C
      if(ier.ne.0) return
C
C
      JTYP = itarg  ! -= zone boundary or zone ctr
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
            Call copyr4(DatBuf(Iwrk2), DatBuf(Iwrk1), INXMR*NTR)
 
C       MAP TO Zone Boundary GRID for each time.
 
      DO 180 IT=1,NTR
 
        IPFL=Iwrk1+(IT-1)*INXMR-1
 
        IPWL=Iwrk2+(IT-1)*INX-1
 
        do ix=1,inx
           if(itarg.eq.2) then
              zinside=DatBuf(IPFL+INX+1-IX)
              zoutsid=DatBuf(IPFL+IX+INX+1)
           else
              zinside=
     >           0.5*(DatBuf(IPFL+INX+1-IX)+DatBuf(IPFL+INX+1-IX+1))
              zoutsid=
     >           0.5*(DatBuf(IPFL+IX+INX+1)+DatBuf(IPFL+IX+INX))
           endif
 
           If (Imap .Eq. 1 ) Then
              DATBUF(IPWL+IX)=zoutsid   ! outside center
 
           Else If (Imap .Eq. 2 ) Then
              DATBUF(IPWL+IX)=zinside   ! Inside center
 
           Else If (Imap .Eq. 3 ) Then
              DATBUF(IPWL+IX)=0.5*(zinside+zoutsid) ! Average in/out
 
           End If
 
        enddo
 
C        LOOP END
 180  CONTINUE
C
C  HAVING REMAPPED DATA, PASS Zone Boundary AS X AXIS TO Calculator
C
      Jtype = JTYP
C
      Return
C
      END
C
C******************** END FILE PLCBND.FOR ; GROUP RPLOT_SUB  ******************
 
