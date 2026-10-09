      Subroutine Findint(Zrecord, Iptr, Int0, Ier)
      use nltrdat_mod

      Character*(*) Zrecord
      Integer        Iptr
      Integer        Ier
      Integer        Int0

      Character*10 Cnumbers

      Data Cnumbers  /'0123456789'/
 
      Int0 = 0    ! Initialize.
      Ier = 0
 
      do 100 I=Iptr,NrecLen
          Ind = Index(Cnumbers, Zrecord(I:I))
          If (Ind .Gt. 0) Then
      	Int0 = Int0*10 + Ind-1
          Else
              Iptr = I+1
      	If (Zrecord(I:I) .Eq. ' ') Then
      	    Go TO 101    ! End Int
      	Else If ( Zrecord(I:I) .Eq. ',') Then
      	    Go To 101    ! End
      	Else
      	    Ier = 1
      	    Go to 101
      	End IF   ! Zrecord
          End If    ! Ind
 
  100   Continue
      Return    ! End of Record
 
  101 Continue
      Return
      End
