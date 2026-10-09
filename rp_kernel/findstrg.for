      Subroutine FindStrg(Zrecord, Iptr, Cstring, Ier)
      use nltrdat_mod
 
      Character*(*) Zrecord
      Character*100  Cin
      Integer        Iptr
      Integer        Ier
      Character*(*)  Cstring
      Character*1    Quote

      Character*10 Cnumbers
      data Cnumbers /'0123456789'/
      data Quote /''''/    ! "'"
 
C	-----------------
 
 
      LenC = LEN(Cstring)
      LenZ = LEN(Zrecord)
      Cstring = ' '  ! Initialize.
      Ier = 0
 
      If (Zrecord(Iptr:Iptr) .Ne. Quote) Then  ! Expect starting quote.
          Ier = 1
          Go to 9000
      End If   ! quote
 
      Iptr = Iptr+1
      Ia   = Iptr-1   ! Current position in input  string.
      Ic = 0          ! Current position of return string.
 
      Do 100 I=Iptr,LenC-1
C	    .Search for ending quote.
            Ia = Ia+1
          If ( LenZ .Lt. Ia) Then
      	Ier = 1
      	Go To 9000
            End If  ! LenZ
          If ( Zrecord(Ia:Ia) .Eq. Quote) Then
C		.See if this is first of a pair or an ending quote.
      	If ( LenZ .Eq. Ia) Then
      	    Iptr = Ia+1
      	    Go to 101
      	End If   ! Lenz
 
      	If ( Zrecord(Ia+1:Ia+1) .Eq. Quote ) Then
C		    .Reduce two quotes to one in output string.
      	    Ia = Ia+1
      	Else
C		    .Ending quote -
      	    Iptr = Ia+1
      	    Go To 101
              End IF    ! Quote
          End if
 
          Ic = Ic+1
          Cstring(Ic:Ic) = Zrecord(Ia:Ia)
 
  100   Continue
 
      If (Zrecord(LenC:LenC) .Ne. Quote) Then
C          .Error - Input line must end with "'" since we have unclosed string.
         Ier = 1
         Go TO 9000
      End If
      Iptr = LenC
 
  101 Continue
 
      Return
 
C	--------------------------
 
 9000   Continue    ! Error return
      Ier = 1
 
      End
