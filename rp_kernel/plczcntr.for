C		                             PLCZCNTR.for in RPLOT_SUB
      SUBROUTINE PLCZCNTR ( STACK, TEMPSTCK, INX, ITYPE)
 
C                  PLCZCNTR takes values in STACK and using the workspace
C                            TEMPSTCK, changes INX zone boundary values to
C			     zone centered values.  Sets ITYPE = 1
C		             which says that the new values in STACK are now
C                            zone centered.
C			     Called from PLCFNXCT.for.
 
      DOUBLE PRECISION STACK(INX)
      REAL TEMPSTCK(INX,2)
 
      INTEGER INX
      INTEGER ITYPE
 
C	------------------------------------------------------------------
 
      DO 200 IX=1,INX
          TEMPSTCK(IX,1) = STACK(IX)            ! Dbl to single precision.
 200  CONTINUE
 
      CALL XINTZC (TEMPSTCK(1,1), TEMPSTCK(1,2), INX)
 
      DO 300 IX=1,INX
          STACK(IX) = TEMPSTCK(IX,2)
 300  CONTINUE
 
      ITYPE = 1        ! zone centered.
 
      RETURN
      END
