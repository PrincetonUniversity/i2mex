      Subroutine PlcSqrt (A)
 
C	.Created 7/26/93 by Tbt. Called from PlcFnxct in RPLOT to get
C	.around a subscript problem which occurs in the handler when
C	.an arithmetic exception happens.
 
C	.Take sqrt of A
      use rpcalc_mod

      Double Precision A

      If ( A .GE. 0.0) Then
          A = Sqrt(A)
      Else
C	    .Error -  - pretend to be handler. See PlcHndlr
          Itype = 6  ! Sqrt of negative number.
          Call PLC_Error (Itype)
          A = 0.0
      End If    ! A<0
	
	
      Return
      End
