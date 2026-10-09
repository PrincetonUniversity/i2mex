C******************** START FILE DATCPU.FOR ; GROUP DATCPU ******************
      SUBROUTINE DATCPU(ZTIME)
C
C  THIS ROUTINE RETURNS THE AMOUNT OF CPU TIME USED UP UNTIL
C  THIS TIME IN SECONDS.
C
C  THIS VERSION FOR USE OUTSIDE OF TRANSP
C
C--------------------------------------
C
      call cpu_time(ztime)
C
      RETURN
      END
C******************** END FILE DATCPU.FOR ; GROUP DATCPU ******************
