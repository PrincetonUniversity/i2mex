      Subroutine FindReal(Zrecord, Iptr, Real0, Ier)
      use nltrdat_mod

      Character*(*) Zrecord
      Character*100  Cin
      Integer        Iptr
      Integer        Ier
      Real           Real0

      Character*10 Cnumbers
      data Cnumbers /'0123456789'/

      Real0 = 0.0    ! Initialize.
      Ier = 0
 
      Do 100 I=Iptr,NrecLen
C	    .Search for delimiter - either " " or ","
          If (Zrecord(I:I) .Eq. ' '  .Or.
     1          Zrecord(I:I) .Eq. ','       ) Then
      	  Ilast = I
      	  Go To 101
          End if
  100   Continue
      Ilast = NrecLen
 
  101 Continue
      	  Cin = Zrecord(Iptr:Ilast)
                Iptr = Ilast
      	  Read ( Cin, 9333, err=8001) Real0
 9333             Format( BN, E15.0)
      Return
 
 
 8001   Continue    ! Error return
      Ier = 1
 
      End
