C                                              PLCHNDLR.FOR IN GROUP RPLOT_SUB
 
      Subroutine PLC_hand
 
 
C	LAST CHANGED:
C	 8/16/93 TBT Created from PLChndlr.for
C	-----------------------------------------------
C
C
C	INTEGER FUNCTION PLCHNDLR IS IN FILE SOURCE:PLCHNDLR.FOR  TBT  7/90
C                                 HANDLES ARITHMETIC EXCEPTIONS FOR RPLOT
C                                 CALCULATOR. STORES 0.0 IS ANY RESULTANT THAT
C                                 HAS AN ARITHMETIC EXCEPTION.
C         Entry Point PLCerror     is called by routines PLCdiv,mult,add,sub
C                                   to simulate error handling & record keeping.
C	  ENTRY POINT PLC_hand IS SUPPLIED TO INITIALIZE # OF ERRORS TO 0.
C	  ENTRY POINT PLC_SUMMARY IS SUPPLIED TO PRINT OUT #'S ERROR MESSAGES.
C	  SUBROUTINE  PLCFIXUP IS CALLED TO ZERO OUT ILLEGAL VALUES PUT INTO THE
C                              RESULTANT AFTER AN ARITHMETIC EXCEPTION OCCURS.
 
C            	         PLCHNDLR IS ESTABLISHED IN PLCFXT.FOR - THE CALCULATOR.
 
      use datmgr_mod
      use rpcalc_mod
 
 
C	-------------------------------------------------------
C	----------------------------------------------  ENTRY POINT PLC_hand
 
C	    .RESET THE # OF ERRORS ENCOUNTERED TO ZERO.
      DO I=1,IAEXS
         INUMERR(I) = 0
      END DO                            ! I
 
      ILERROR = .FALSE.                 ! NO ERROR YET FOR NEW CACULATION
      IERRPOS = 1
      INDEX1  = 0
 
      RETURN
 
      END
