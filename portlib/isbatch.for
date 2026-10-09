      subroutine isbatch(iret)
C
C  DMC 9 SEP 1991
C
c
c  this fortran routine replaces the old macro "isbatch".  This routine
c  returns IRET=1 if the current process is a batch job or spawned
c  subprocess not receiving input from a terminal; it returns IRET=0
c  if an interactive terminal is in control.
c
      INTEGER FISATTY
C
      IANS=FISATTY(0)
C
      IF(IANS.EQ.1) THEN
        IRET = 0
      ELSE
        IRET = 1
      ENDIF
C
 9000 continue
      return
      end
