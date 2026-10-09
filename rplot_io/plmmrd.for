C-----------------------------------------------------------------------
C  PLMMRD -- READ MOMENTS DATA
 
 
      SUBROUTINE PLMMRD(IDIM,IBOOST)
 
C  dmc 4 Aug 1999 -- changed UERMSG calls to ZERMSG calls -- avoid
C  UREADSUB dependence.
C    These error messages won't occur for existing TRANSP runs if
C  the "MOMRUN" test passes.
C
C---------------------
C  PASSED CONTROLS
C   IDIM = 1 TO READ LIMITED INFO ONLY (FOR MP LINE INTEGRAL,
C             CF PLMMDR SUBROUTINE)
C        = 2 TO READ R,Y MOMENTS -- COMPLETE DATA SET
C
C   IBOOST = 1 TO LEAVE STORAGE PRIORITY OF ALL MOMENTS DATA BOOSTED TO
C            "6" ON EXIT; = 0 TO RESTORE STANDARD PRIO OF "5" ON EXIT
C
C  OUTPUT
C    PLFMPA COMMON -- PTRS TO FUNCTIONS READ IN
C
C  READ DATA TO MEMORY BUFFER; OUTPUT ADDRESSES TO PLFMPA COMMON
 
 
C 04/13/94 tbt Assigned dummy functions to zero moment sin terms.
 
 
 
C
C  COMMON BLOCKS---
C
      use datmgr_mod
      use cplotr_mod
      use plfmpa_mod
C
C  LOCAL--
C
      LOGICAL MOMRUN,FUNGOT
C
C-----------------------------------------------------------------------
C
      IMOM=0
C
      IF(.NOT.MOMRUN(INUM,IMMGEO)) THEN
        CALL ZERMSG(
     >'? PLOTMM -- MOMENTS EQUILIBRIUM DATA NOT AVAILABLE, '//RUNLB2)
        RETURN
      ENDIF
C
C  READ THE 0TH R MOMENT DATA
C   0TH R AND Y MOMENTS FOR ASYMMETRIC CASES
C
      IF(IMMGEO.EQ.0) THEN
        CALL DMGFXT(INUM,INDR0)
        CALL PLPRIO(6,INDR0)
        ILR0=INUM
      ENDIF
C
C  Asymmetric:
C
      IF((IMMGEO.EQ.1).AND.(IDIM.EQ.2)) THEN
        PLTABB='RMC00'
        IF(.not.FUNGOT(INUM)) then
          CALL ZERMSG('?PLMMRD:  MISSING REQUIRED FCN:  '//PLTABB)
        ENDIF
        CALL DMGFXT(INUM,INDF)
        CALL PLPRIO(6,INDF)
        INDAMM(1,0)=INDF
        ILCAMM(1,0)=INUM
C
        ILR0=INUM  ! for labeling...
C
        INDAMM(2,0)=Indf
        ILCAMM(2,0)=Inum            ! Dummy assign  tbt
C
        PLTABB='YMC00'
        IF(.not.FUNGOT(INUM)) then
          CALL ZERMSG('?PLMMRD:  MISSING REQUIRED FCN:  '//PLTABB)
        ENDIF
        CALL DMGFXT(INUM,INDF)
        CALL PLPRIO(6,INDF)
        INDAMM(3,0)=INDF
        ILCAMM(3,0)=INUM
C
        INDAMM(4,0)= Indf           ! Dummy assign    tbt
        ILCAMM(4,0)= Inum           ! Always multiplied by sin(0.)
      ENDIF
C
      IF(IMMGEO.EQ.1) THEN
        PLTABB='YMPA'
        IF(.not.FUNGOT(INUM)) then
          CALL ZERMSG('?PLMMRD:  MISSING REQUIRED FCN:  '//PLTABB)
        ENDIF
        CALL DMGFXT(INUM,INDF)
        CALL PLPRIO(6,INDF)
        INDYMP=INDF
        ILYMP=INUM
C
        PLTABB='RMAJM'
        IF(.not.FUNGOT(INUM)) then
          CALL ZERMSG('?PLMMRD:  MISSING REQUIRED FCN:  '//PLTABB)
        ENDIF
        CALL DMGFXT(INUM,INDF)
        CALL PLPRIO(6,INDF)
        INDRMP=INDF
        ILRMP=INUM
      ENDIF
C
C  READ HIGHER ORDER MOMENTS
C   MOMENTS ARE ASSUMED TO EXIST IN PAIRS; I.E. IF THERE IS A 3RD R
C  MOMENT THEN THERE IS ALSO A 3RD Y MOMENT
C
      IMOM=1
      IF((IMMGEO.EQ.1).AND.(IDIM.EQ.1)) GO TO 20
C
      IF(IMMGEO.EQ.1) THEN
C-------------------------------------------------
C  ASYMMETRIC MOMENTS SET
C  mod dmc -- deal with variation in moments profile names
C    RMC010 or RMC10 would both be names for the cos(10*theta) R moment...
C
        DO IMOM=1,NAXMOM
          WRITE(PLTABB,'(''RMC'',I2.2)') IMOM
          IF(.NOT.FUNGOT(IFNRC)) then
             WRITE(PLTABB,'(''RMC0'',I2.2)') IMOM
             IF(.NOT.FUNGOT(IFNRC)) GO TO 20
          endif
          WRITE(PLTABB,'(''RMS'',I2.2)') IMOM
          IF(.NOT.FUNGOT(IFNRS)) then
             WRITE(PLTABB,'(''RMS0'',I2.2)') IMOM
             IF(.NOT.FUNGOT(IFNRS)) GO TO 20
          endif
          WRITE(PLTABB,'(''YMC'',I2.2)') IMOM
          IF(.NOT.FUNGOT(IFNYC)) then
             WRITE(PLTABB,'(''YMC0'',I2.2)') IMOM
             IF(.NOT.FUNGOT(IFNYC)) GO TO 20
          endif
          WRITE(PLTABB,'(''YMS'',I2.2)') IMOM
          IF(.NOT.FUNGOT(IFNYS)) then
             WRITE(PLTABB,'(''YMS0'',I2.2)') IMOM
             IF(.NOT.FUNGOT(IFNYS)) GO TO 20
          endif
C
          CALL DMGFXT(IFNRC,IND)
          CALL PLPRIO(6,IND)
          INDAMM(1,IMOM)=IND
          ILCAMM(1,IMOM)=IFNRC
C
          CALL DMGFXT(IFNRS,IND)
          CALL PLPRIO(6,IND)
          INDAMM(2,IMOM)=IND
          ILCAMM(2,IMOM)=IFNRS
C
          CALL DMGFXT(IFNYC,IND)
          CALL PLPRIO(6,IND)
          INDAMM(3,IMOM)=IND
          ILCAMM(3,IMOM)=IFNYC
C
          CALL DMGFXT(IFNYS,IND)
          CALL PLPRIO(6,IND)
          INDAMM(4,IMOM)=IND
          ILCAMM(4,IMOM)=IFNYS
        ENDDO
        IMOM=NAXMOM+1
        GO TO 20
C-------------------------------------------------
      ENDIF
C-------------------------------------------------
C  SYMMETRIC MOMENTS SET
      DO 10 IMOM=1,NAXMOM
C
        PLTABB=' '
        WRITE(PLTABB,1001) IMOM
 1001   FORMAT('RMM',I2.2)
        IF(.NOT.FUNGOT(IFNR)) then
           WRITE(PLTABB,'(''RMM0'',I2.2)') IMOM
           IF(.NOT.FUNGOT(IFNR)) GO TO 20
        endif
C
        IF(IDIM.GT.1) THEN
          WRITE(PLTABB,1002) IMOM
 1002     FORMAT('YMM',I2.2)
          IF(.NOT.FUNGOT(IFNY)) then
             WRITE(PLTABB,'(''YMM0'',I2.2)') IMOM
             IF(.NOT.FUNGOT(IFNY)) GO TO 20
          endif
        ENDIF
C
C  ANOTHER MOMENT EXISTS; READ R THEN Y
        CALL DMGFXT(IFNR,IND)
        CALL PLPRIO(6,IND)
        INDMM(1,IMOM)=IND
        ILCMM(1,IMOM)=IFNR
C
        IF(IDIM.GT.1) THEN
          CALL DMGFXT(IFNY,IND)
          CALL PLPRIO(6,IND)
          INDMM(2,IMOM)=IND
          ILCMM(2,IMOM)=IFNY
        ENDIF
C
C  END LOOP
 10   CONTINUE
C
      IMOM=NAXMOM+1
C-------------------------------------------------
 20   CONTINUE
C
C  NUMBER OF AVAILABLE MOMENTS
      IMOM=IMOM-1
C
C  RESTORE NORMAL PRIORITIES UNLESS BOOST FLAG IS SET
C
      IF(IBOOST.NE.1) THEN
C
        CALL PLMMRLS(IDIM)
C
      ENDIF
C
      RETURN
      END
