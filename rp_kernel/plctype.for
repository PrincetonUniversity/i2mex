C						PLCTYPE.for in RPLOT_SUB
 
      SUBROUTINE PLCTYPE (STACK1, STACK2, TEMPSTCK, INX,
     A                      ITYPE,  ZINPUT, IT, ICURPOS, LERR, LSHIFT)
 
C                  PLCTYPE checks to make sure that the TYPEs of the x-axis
C                          of the operands are the same. If one is zone centered
C                          and the other is zone boundary, the zone boundary one
C                          is transformed to zone centered.
C 		           Called by PLCFNXCT.
C
C  Modified:
C  06/01/95 TBT Comment out call to uermsg.
 
 
C	Arguments input
C			     STACK1 & STACK2 - input stack values for 2 operands
C			     TEMPSTCK - workspace array.
C			     INX - # of values for the X-axis.
C			     ITYPE - I array with X-axis type for 2 operands.
C			     ZINPUT - C*128 input string to parser. Used for
C                                     error messages.
C			     IT - I*4 ITth time. Print errors for IT=1 only.
C                            ICURPOS - I*4 current cursar position in input line
C	Argument returned.
C			     STACK1 & STACK2 can change to zone centered.
C	                     LERR - Logical set to true if major X-axis
C  			             inconsistency.
C                            LSHIFT set to .TRUE. if zone-ctr shift is
C                                    done.  LSHIFT is not cleared here.
C
 
 
      DOUBLE PRECISION STACK1(INX), STACK2(INX)
      REAL   TEMPSTCK(INX,2)
      INTEGER ITYPE(*)
      LOGICAL LERR,LSHIFT
      LOGICAL LXOUT
      CHARACTER*(*) ZINPUT
 
      logical transp_run_ck
 
      CHARACTER*1 HYPHEN
      CHARACTER*1 CARROT
 
      data hyphen/'-'/
      data carrot/'^'/
 
C  	---------------------------------------------------------------------
 
 
      LERR = .FALSE.
 
      IF (ITYPE(1) .LE. 0) THEN
         ITYPE(1) = ITYPE(2)	    ! First operand is a constant or scalar.
         GO TO 9000
      END IF
 
      IF (ITYPE(2) .LE. 0)  GO TO 9000  ! Second operand is a cnstnt or sclr.
 
      IF (ITYPE(1) .EQ. ITYPE(2)) GO TO 9000  ! X-axis agree
 
      if(transp_run_ck(0)) then
         IF (ITYPE(1) .EQ. 1  .AND.  ITYPE(2) .EQ. 2)  THEN
C	    .Second operand is zone boundary and first is zone centered.
C	    .Change second to zone centered.
	
            LSHIFT = .TRUE.
            CALL PLCZCNTR( STACK2, TEMPSTCK, INX, ITYPE(2))
            IF (IT .EQ. 1)  GO TO 8000  ! Print warning for time = 1.
            GO TO 9000
         END IF                         ! 1&2
 
         IF (ITYPE(1) .EQ. 2  .AND. ITYPE(2) .EQ. 1)  THEN
C	    .First operand is zone boundary and second is zone centered.
C	    .Change first to zone centered.
	
            LSHIFT = .TRUE.
            CALL PLCZCNTR( STACK1, TEMPSTCK, INX, ITYPE(1))
            IF (IT .EQ. 1)  GO TO 8000  ! Print warning for time = 1.
            GO TO 9000
         END IF                         ! 2&1
      endif
 
 
      LERR = .TRUE.
      GO TO 9000
 
C	.Print warning message for all zone changes for time = 1.
 8000 continue
 
 9000 CONTINUE
 
 
      RETURN
      END
 
 
      logical function transp_run_ck(idum)
      use cplotr_mod
 
      transp_run_ck = NLTRANSP
 
      return
      end
