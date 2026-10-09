C******************** START FILE TFILW2.FOR ; GROUP TFILIO ******************
C---------------------------------------------------------------
C  TFILW2
C
C  WRITE TRANSP PLOTTING LABEL INFORMATION
C
C  SUBROUTINE OF TFILWR
C
C  BROKEN OUT DMC 16 JUNE 1988
C  mod DMC December 2009 -- F09 format supported
C    to allow 64 char labels and 32 units labels
C    at present abbreviations are still 10 characters but could be expanded 
C    to 32 with minimal adjustment to tfilw2 & tfilr2
C      (A,22X) -> (A,12X) for 20 char abbrev; (A,22X) -> (A) for 32 chars.
C
C------------------------------
C
      SUBROUTINE TFILW2(LUN,LABEL,UNITS,ABBREV,IARG1,IARG2,ITYPE)
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
C  ITYPE=1 -- SCALAR FUNCTION
C  ITYPE=2 -- PROFILE FUNCTION
C  ITYPE=3 -- MULTIGRAPH
C
C---------------------------------------------------------------
C
C  EXECUTABLE CODE
C
 2099 format(6x,a)
C
      IF(ITYPE.EQ.1) THEN
C  SCALAR FUNCTION LABEL
         write(lun,2011) abbrev,iarg1
         write(lun,2099) units
         write(lun,2099) label
 2011    FORMAT(A,t35,i5)
C
      ELSE IF(ITYPE.EQ.2) THEN
C  PROFILE FUNCTION LABEL
         IF(NLXVAR) THEN
            write(lun,2021) abbrev,iarg1,iarg2
            write(lun,2099) units
            write(lun,2099) label
 2021       FORMAT(A,t35,i5,1x,i5)
         ELSE
            WRITE(LUN,2022) ABBREV,IARG1
            write(lun,2099) units
            write(lun,2099) label
 2022       FORMAT(A,t35,i5)
         ENDIF
C
      ELSE IF(ITYPE.EQ.3) THEN
C  MULTIGRAPH LABEL
         WRITE(LUN,2031) abbrev,iarg1,iarg2
         write(lun,2099) units
         write(lun,2099) label
 2031    FORMAT(A,t35,i5,1x,i5)
      ENDIF
C
      RETURN
      END
C******************** END FILE TFILW2.FOR ; GROUP TFILIO ******************
