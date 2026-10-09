      Subroutine PlcMult(A, B)
 
C	.Created 11/6/92 by Tbt. Called from PlcFnxct in RPLOT to get
C	.around a subscript problem which occurs in the handler when
C	.an arithmetic exception happens.
 
C	.Multiply two numbers.
 
C	Updates:
C	07/26/93 TBT Added checks for overflow and call to PLC_error.
      use rpcalc_mod

      Double Precision A,B

      Real ZlogA, ZlogB
      Integer  Itype
 
 
      If ( A .Eq. 0.0D0 .OR. B .Eq. 0.0D0) Then
         A = 0.0D0
      Else
         ZlogA = Log10( Abs(A))
         ZlogB = Log10( Abs(B))
 
         If ((ZlogA+ZlogB) .Lt. ExpMax ) Then
            If ((ZlogA+ZlogB) .Gt. ExpMin ) Then
               A = A * B
            Else
C     .Probable underflow  -- silently set to zero
csilent               Itype = 3
csilent               Call PLC_error (Itype)
               A = 0.0
            End If                      ! Expmin
         Else
C		.Probable overflow.
            Itype = 1
            Call PLC_error (Itype)
            A = 0.0
         End If                         ! Expmax
      End If                            ! 0.0
	
      Return
      End
