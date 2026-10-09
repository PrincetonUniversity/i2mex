      Subroutine PlcExp (A)
 
C	.Created 11/6/92 by Tbt. Called from PlcFnxct in RPLOT to get
C	.around a subscript problem which occurs in the handler when
C	.an arithmetic exception happens.
 
C	.Find EXP(A)
      use rpcalc_mod

      Double Precision A

      If ( A .LT. (Expmax / log10(exp(1.d0))) ) Then
          A = EXP(A)
      Else
C	    .Error - A too big  - pretend to be handler. See PlcHndlr
          Itype = 1  ! Divide by zero.
          Call PLC_Error (Itype)
          A = 0.0D0
      End If
	
	
      Return
      End
