C                                              PLC_Error.FOR IN GROUP RPLOT_SUB
 
      Subroutine PLC_error (Itype)
 
C       _______ ________ ________
C
C	LAST CHANGED:
C        8/18/93 TBT Set itype = itype (was = 1)
C	 8/16/93 TBT Created from PLChndlr.for
C	-----------------------------------------------
C
 
C	INTEGER FUNCTION PLCHNDLR IS IN FILE SOURCE:PLCHNDLR.FOR  TBT  7/90
C                                 HANDLES ARITHMETIC EXCEPTIONS FOR RPLOT
C                                 CALCULATOR. STORES 0.0 IS ANY RESULTANT THAT
C                                 HAS AN ARITHMETIC EXCEPTION.
C         Entry Point PLCerror     is called by routines PLCdiv,mult,add,sub
C                                   to simulate error handling & record keeping.
C	  ENTRY POINT PLC_HAND_INI IS SUPPLIED TO INITIALIZE # OF ERRORS TO 0.
C	  ENTRY POINT PLC_SUMMARY IS SUPPLIED TO PRINT OUT #'S ERROR MESSAGES.
C	  SUBROUTINE  PLCFIXUP IS CALLED TO ZERO OUT ILLEGAL VALUES PUT INTO THE
C                              RESULTANT AFTER AN ARITHMETIC EXCEPTION OCCURS.
 
C            	         PLCHNDLR IS ESTABLISHED IN PLCFXT.FOR - THE CALCULATOR.
 
	
C
      use datmgr_mod
      use rpcalc_mod
 
      Integer Itype
 
      Logical LXOUT
 
C	----------------------------------------------  Entry point PLC_error
C	.Do accounting for a found error in another routine.
 
 
C	    Input argument - Itype - Integer giving index of Carthex type
C                                    of error.
 
      If (Itype .Le. 0  .Or. Itype .Gt. Iaexs) Then
         Write (lunzer(0), 9000) Itype
 9000    Format(/' Plc_error: Itype out of range =', I5/)
         call abortt
      End If                            ! Lxout
 
C	        Abort !!!
 
	
      Inumerr(Itype) = Inumerr(Itype) + 1
 
      IF (.NOT. ILERROR) THEN           !
C		. THIS IS THE FIRST ERROR FOR THIS INPUT LINE.	!
         ILERROR = .TRUE.               !
         IERRPOS = ICURPOS              ! CURRENT POSITION IS FIRST ERROR POSITION
         INDEX1  = Itype                ! INDEX OF FIRST ERROR.
      END IF                            ! ILERROR
 
      Return
 
      END
