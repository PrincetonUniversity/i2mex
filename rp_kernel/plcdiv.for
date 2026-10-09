      Subroutine PlcDiv (A, B)
 
C	.Created 11/6/92 by Tbt. Called from PlcFnxct in RPLOT to get
C	.around a subscript problem which occurs in the handler when
C	.an arithmetic exception happens.
 
C	Updates:
C	 8/16/93  TBT  Added underflow.
 
 
C	.Divide two numbers.
      use rpcalc_mod

      Double Precision A,B
      Double Precision ZLogA, ZlogB

      If ( B .Ne. 0.0) Then
         If ( A .Eq. 0.0) Then
            A = 0.0
         Else                           ! A & B are nonzero.
            ZLogA = Log10(Abs(A))
            ZLogB = Log10(Abs(B))
            If (ZLogA-ZLogB .Lt. ExpMax) Then
               If (ZlogA-ZlogB .Gt. Expmin) Then
                  A = A / B
               Else
C     .Error - floating point underflow.  -- silently set to zero
csilent                  Itype = 3
csilent                  Call PLC_error (Itype)
                  A = 0.0
               End If                   ! Underflow
            Else
C     .Error - floating pt overflow.
               Itype = 1
               Call PLC_error (Itype)
               A = 0.0
            End If                      ! ZlogA
         End If                         ! A=0
 
      Else
C	    .Error - divide by zero - pretend to be handler. See PlcHndlr
         Itype = 2                      ! Divide by zero.
         Call PLC_Error (Itype)
         A = 0.0
      End If                            ! b=0
	
	
      Return
      End
