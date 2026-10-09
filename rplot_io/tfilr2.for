C******************** START FILE TFILR2.FOR ; GROUP TFILIO ******************
C---------------------------------------------------------------
C  TFILR2
C
C  READ TRANSP PLOTTING LABEL INFORMATION
C
C  SUBROUTINE OF TFILRD
C
C  BROKEN OUT DMC 16 JUNE 1988
C  mod DMC December 2009 -- F09 format supported
C    to allow 64 char labels and 32 units labels
C    at present abbreviations are still 10 characters but could be expanded 
C    to 32 with minimal adjustment to tfilw2 & tfilr2
C      (A,22X) -> (A,12X) for 20 char abbrev; (A,22X) -> (A) for 32 chars.
C
      SUBROUTINE TFILR2(LUN,LABEL,UNITS,ABBREV,IARG1,IARG2,ITYPE,FMT)
C
C  COMMON BLOCKS--
C
      use cplotr_mod
C
C  PASSED--
C
      INTEGER LUN  ! FILE LUN, FILE IS ALREADY OPEN!
C
      CHARACTER*(*) LABEL    ! FUNCTION OR MULTIGRAPH LABEL
      CHARACTER*(*) UNITS    ! FUNCTION OR MULTIGRAPH UNITS
      CHARACTER*(*) ABBREV   ! FUNCTION OR MULTIGRAPH ABBREVIATION
C
      INTEGER IARG1          ! AUXILLIARY INTEGER 1
      INTEGER IARG2          ! AUXILLIARY INTEGER 2
C
      INTEGER ITYPE          ! TYPE OF ITEM BEING LABELED
C
      CHARACTER*3 FMT        ! FORMAT CODE, '   ' = OLD STYLE
C
C  ITYPE=1 -- SCALAR FUNCTION
C  ITYPE=2 -- PROFILE FUNCTION
C  ITYPE=3 -- MULTIGRAPH
C
C  LOCAL
C
      CHARACTER*1 ZCT
C
C---------------------------------------------------------------
C
C  EXECUTABLE CODE
C
      LABEL=' '
      UNITS=' '
      ABBREV=' '
      iarg1=0
      iarg2=0
C
 2099 format(6x,a)
C
      IF(ITYPE.EQ.1) THEN
C  SCALAR FUNCTION LABEL
        IF(FMT.EQ.'   ') THEN
          READ(LUN,2011) LABEL(1:20),UNITS(1:10),ABBREV(1:5),IARG1
        ELSE IF(FMT.EQ.'F88') THEN
          READ(LUN,2011) LABEL(1:32),UNITS(1:16),ABBREV,IARG1
        ELSE IF(FMT.EQ.'F09') THEN
          READ(LUN,'(A,t35,i5)') abbrev,iarg1
          read(lun,2099) units
          read(lun,2099) label
        ENDIF
 2011   FORMAT(A,A,A,I1)
C
      ELSE IF(ITYPE.EQ.2) THEN
C  PROFILE FUNCTION LABEL
        IF(NLXVAR) THEN
          IF(FMT.EQ.'   ') THEN
            READ(LUN,2021) LABEL(1:20),UNITS(1:10),ABBREV(1:5),
     >          IARG1,ZCT,IARG2
          ELSE IF(FMT.EQ.'F88') THEN
            READ(LUN,2021) LABEL(1:32),UNITS(1:16),ABBREV,IARG1,ZCT,
     >            IARG2
          ELSE IF(FMT.EQ.'F09') THEN
            ZCT=' '
            READ(LUN,'(A,t35,i5,1x,i5)') abbrev,iarg1,iarg2
            read(lun,2099) units
            read(lun,2099) label
          ENDIF
C  KLUGE FOR I1 VS I2 FMT IN DIFFERENT REVS OF LABEL FILE
 2021     FORMAT(A,A,A,I1,A1,I4)
          IF(ZCT.NE.' ') THEN
             READ(ZCT,'(I1)') ITEMP
             IARG1=10*IARG1+ITEMP
          ENDIF
        ELSE
C  PROFILE FCN LABEL, FIXED X AXIS FORMAT
          IF(FMT.EQ.'   ') THEN
            READ(LUN,2022) LABEL(1:20),UNITS(1:10),ABBREV(1:5),IARG1
          ELSE IF(FMT.EQ.'F88') THEN
            READ(LUN,2022) LABEL(1:32),UNITS(1:16),ABBREV,IARG1
          else if(FMT.EQ.'F09') THEN
            read(lun,'(A,t35,i5)') abbrev,iarg1
            read(lun,2099) units
            read(lun,2099) label
          ENDIF
 2022     FORMAT(A,A,A,I1)
        ENDIF
C
      ELSE IF(ITYPE.EQ.3) THEN
C  MULTIGRAPH LABEL
        IF(FMT.EQ.'   ') THEN
          READ(LUN,2031) LABEL(1:20),UNITS(1:10),
     >        IARG1,IARG2,ABBREV(1:5)
        ELSE IF(FMT.EQ.'F88') THEN
          READ(LUN,2031) LABEL(1:32),UNITS(1:16),IARG1,IARG2,ABBREV
        else if(FMT.EQ.'F09') THEN
           read(lun,'(A,t35,i5,1x,i5)') abbrev,iarg1,iarg2
           read(lun,2099) units
           read(lun,2099) label
        ENDIF
 2031   FORMAT(A,A,2I5,A)
      ENDIF
C
C  REMOVE "$" FROM LABELS
C
      CALL TFILND(UNITS)
C
      RETURN
      END
C******************** END FILE TFILR2.FOR ; GROUP TFILIO ******************
