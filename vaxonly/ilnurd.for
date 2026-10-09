C******************** START FILE ILNURD.FOR ; GROUP URDUTS ******************
C--------------------------------------------------------------
C&  ILNURD  INTEGER FCN (UREAD) LEN OF STRING W/O TRAILING BLANKS
C
      INTEGER FUNCTION ILNURD(STR)
C
      implicit NONE
C
      CHARACTER*(*) STR
C
      integer str_length
C
C---------------------------------------
C  replace legacy code with portlib call -- dmc Oct 2002
C
      ilnurd=max(1,str_length(str))
C
      END
C******************** END FILE ILNURD.FOR ; GROUP URDUTS ******************
