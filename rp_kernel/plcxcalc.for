C			                      PLCXCALC.FOR in RPLOT_SUB
 
      SUBROUTINE PLCXCALC( STACK, TEMPSTCK, INX, IINT, IT)
 
C	           PLCXCALC sets up the call to IXCALC to do various integrals.
C		            It is called from PLCFNXCT.    TBT 9/27/90
 
C	Arguments:
C	    Input:
C		STACK(INX) - R*8 array containing values for each point on the
C                            X axis to be evaluated by the function.
C               TEMSTCK(*,2) - R*4 array to be used as a work space. Wiped out.
C	        INX          - I*4 number of points in the X direction.
C               IINT         - I*4 code giving which function to use on STACK
C                                  data. See IXCALC for codes.
C		IT           - I*4 ITth time.
C	    Output:
C		STACK(INX)   - R*8 array of function of input value.
 
 
 
      use cplotr_mod
 
      DOUBLE PRECISION  STACK(INX)
 
      REAL   TEMPSTCK( NR0, 2)    ! TEMPORARY STACK.
 
C	----------------------------------------------------------------
 
      DO 110 IX=1,INX
          TEMPSTCK(IX, 1) = STACK(IX)   ! Double to single precision.
 110  CONTINUE
 
      ZSGN = 1.
      ZX0L = 0
      ZF0  = 0
 
      IINTL = IINT
      CALL DMGGEO( IINTL, IND8, IND9)   ! GET GEOMETRY INFORMATION.
					  !   IND8 & IND9 ARE DUMMIES.
 
      CALL IXCALC( TEMPSTCK, TEMPSTCK(1,2), INX, ZSGN,
     1               IINT, IT, ZX0L, ZF0                )
 
      DO 210 IX=1,INX
          STACK(IX) = TEMPSTCK(IX, 2)   ! Single to double precision.
 210  CONTINUE
 
      RETURN
      END
 
 
