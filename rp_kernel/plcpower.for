      Subroutine PlcPower (A, B)
 
C	.Created 7/26/93 by Tbt. Called from PlcFnxct in RPLOT to get
C	.around a subscript problem which occurs in the handler when
C	.an arithmetic exception happens.
 
C	.Raise A to the B power.
      use rpcalc_mod

      Double Precision A,B
      Double Precision dZero

      dZero = 0.0D1
 
      If (A .Eq. dZero .AND. B .Eq. dZero) Then
C	    .Undefined exponentiation.
          Itype = 4
          Call PLC_Error (Itype)
          A = dZero
 
      Else if (A .Eq. dZero) Then
          A = dZero    ! B nonzero.
 
      Else if (B .Eq. dZero) Then
          A = 1.0d0
 
      Else if( Log10(Abs(A))*B .Gt. Expmax ) Then
C	    .Errorr A**B is too large for machine
          Itype = 1  ! Overflow.
          Call PLC_Error (Itype)
          A = dZero
 
      Else
          If (Mod(B, 1.D0) .Ne. 0.D0) Then
C	      .B is not a whole number
            If (A .Gt. 0.0) Then
      	A = A ** B
            Else
C		.Undefined exponentiation - A negative, B not whole #
      	Itype = 4
      	Call PLC_Error (Itype)
      	A = dZero
            End If   ! Mod B
          Else
              A = A ** Int(B)    ! B is whole number
          End If    ! Mod
      End If    !
	
      Return
      End
