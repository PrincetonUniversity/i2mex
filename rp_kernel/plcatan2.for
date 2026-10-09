      Subroutine PlcATan2(A, B)
 
C	.Created 11/6/92 by Tbt. Called from PlcFnxct in RPLOT to get
C	.around a subscript problem which occurs in the handler when
C	.an arithmetic exception happens.
 
C	.Take ATan2 two numbers.
 
C	Updates:
C	07/26/93 TBT Added checks for overflow and call to PLC_error.
      use rpcalc_mod

      Double Precision A,B

      Integer  Itype
 
 
      If ( A .Eq. 0.0D0  .And.  B .Eq. 0.0D0) Then
C	    .Illegal arguments.
          Itype = 10
          Call PLC_error (Itype)
          A = 0.0D0
      Else
          A = Atan2(A,B)
      End If   ! 0.0
	
      Return
      End
