C******************** START FILE FUNGOT.FOR ; GROUP PLOTMM ******************
C-----------------------------------------------------------------
C
      LOGICAL FUNCTION FUNGOT(INUM)
C
C  THIS FCN RETURNS .TRUE. IF THE TRANSP RUN OUTPUT DATA CONTAINS THE
C  NAMED FUNCTION (VS. T + ADDL COORD)
C  .FALSE. OTHERWISE
C
C  THE NAMED FUNCTION IS SPECIFIED IN COMMON VARIABLE PLTABB
C
C  COMMON BLOCKS---
C
      use datmgr_mod
      use cplotr_mod
C---
      CHARACTER*1 ZBUFF(80)
C
      ILENA=LEN(PLTABB)
C
      ZBUFF(ILENA+1)=' '
      DO IC=1,ILENA
        ZBUFF(IC)=PLTABB(IC:IC)
      ENDDO
C
      FUNGOT=.FALSE.
      INUM=ISCMP0(ZBUFF,ABR,NFXT,ILENA,IDUM)
C
      IF(INUM.GT.0) FUNGOT=.TRUE.
C
      RETURN
      END
C******************** END FILE FUNGOT.FOR ; GROUP PLOTMM ******************
