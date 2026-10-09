      Subroutine Findeq(Zrecord, Iptr, Ier)
      use nltrdat_mod
 
      Character*(*) Zrecord
      Integer        Iptr
      Integer        Ier

      Ier = 0
      LenZ = Len(Zrecord)
 
      Do 100 I=Iptr,LenZ
          If (Zrecord(I:I) .Eq. '=') Then
      	Iptr = I
      	Go To 101
          End If    ! =
  100 Continue
 
      Ier = 5     ! no equal sign
      Iptr = 0
      Return
 
  101 Continue
	
 
      Return
      End
