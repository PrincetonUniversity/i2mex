      SUBROUTINE MCRDR_OPEN(LUN,NAME)
C
C  OPEN A FILE - READONLY
C
      CHARACTER*(*) NAME
C
      ILUN=IABS(LUN)
      CLOSE(UNIT=ILUN,ERR=5)
 5    CONTINUE
C
      IF(LUN.GT.0) THEN
        OPEN(UNIT=LUN,FILE=NAME,
     >         status='OLD',ACCESS='SEQUENTIAL')
      ELSE
C
C NONAME OPEN
C
         OPEN(UNIT=ILUN,
     >         status='OLD',ACCESS='SEQUENTIAL')
      ENDIF
C
      RETURN
      END
 
 
 
 
 
 
 
 
