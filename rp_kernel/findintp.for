      Subroutine FindintP(Zrecord, Iptr, Int0, Ier)
      use nltrdat_mod

C	Find an integer foillowed by a ')'

      Character*(*) Zrecord
      Integer        Iptr
      Integer        Ier
      Integer        Int0
 
 
      Character*10 Cnumbers
      data Cnumbers /'0123456789'/
 
      Int0 = 0    ! Initialize.
      Ier = 0
 
      do 100 I=Iptr,NrecLen
          Ind = Index(Cnumbers, Zrecord(I:I))
          If (Ind .Gt. 0) Then
      	Int0 = Int0*10 + Ind-1
          Else
              Iptr = I+1
      	If (Zrecord(I:I) .Eq. ') ') Then
                  If (Zrecord(I+1:I+1) .Ne. '=') Then
                  Ier=1
                End If
      	  Go TO 101    ! End Int
C		Else If ( Zrecord(I:I) .Eq. ',') Then
C		    Go To 101    ! End
      	Else
      	    Ier = 1
      	    Go to 101
      	End IF   ! Zrecord
          End If    ! Ind
 
  100   Continue
 
      Ier = 1   ! No ')' found.
      Return    ! End of Record
 
  101 Continue
      Return
      End
