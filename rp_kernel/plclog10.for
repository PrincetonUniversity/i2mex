      Subroutine PlcLog10 (A)
 
C	.Created 11/6/92 by Tbt. Called from PlcFnxct in RPLOT to get
C	.around a subscript problem which occurs in the handler when
C	.an arithmetic exception happens.
 
C	.Find LOG10(A)
      use rpcalc_mod

      Double Precision A

      If ( A .GT. 0.0d1) Then
          A =Log10(A)
      Else
C	    .Error -  - pretend to be handler. See PlcHndlr
          Itype = 5  ! Log of 0 or negative number.
          Call PLC_Error (Itype)
          A = 0.0
      End If
	
	
      Return
      End
