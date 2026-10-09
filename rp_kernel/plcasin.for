      Subroutine PlcASin (A)
 
C	.Created 7/27/93 by Tbt. Called from PlcFnxct in RPLOT to get
C	.around using a handler ( - make machine independent!!) when
C	.an arithmetic exception happens.
 
C	.Find ASin(A)
      use rpcalc_mod

      Double Precision A

      If ( Abs(A) .LE. 1.0d0) Then
          A = Asin(A)
      Else
C	    .Error -  pretend to be handler. See PlcHndlr
          Itype = 10
          Call PLC_Error (Itype)
          A = 0.0
      End If
	
	
      Return
      End
