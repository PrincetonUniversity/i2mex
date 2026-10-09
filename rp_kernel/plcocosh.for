      Subroutine PlcOCosH (A)
 
C	.Created 11/6/92 by Tbt. Called from PlcFnxct in RPLOT to get
C	.around a subscript problem which occurs in the handler when
C	.an arithmetic exception happens.
 
C	.Find 1/CosH(A)
      use rpcalc_mod

      Double Precision A

      If ( Abs(A) .LT. (Expmax / log10(exp(1.d0))) ) Then
          A = 1.d0/CosH(A)
      Else
C	    .Error - A too big  - pretend to be handler. See PlcHndlr
          Itype = 10  ! Bad argument in math lib.
          Call PLC_Error (Itype)
          A = 0.0D0
      End If
	
	
      Return
      End
