      Subroutine Findlog(Zrecord, Iptr, Nlogical, Ier)
      use nltrdat_mod

C	Find an logical followed by a ','.
C	Return the logical in NLogical
C	A logical is ".T...", "T...", ".F...", or "F..." where
C                      "..." is anything except "," or " " which are delimiters.
 
      Character*(*) Zrecord
      Integer        Iptr
      Integer        Ier
      Logical	       NLogical
 
      Character*1    ThisChar

C	-----------------------
 
      IzLen = Len(Zrecord)
      Ier = 0
      Iq = Iptr-1
 
        Do 500 Ia=Iptr,IzLen
 
          Iq = Iq+1
          If (Iq .Gt. Izlen) Go to 500
          ThisChar = Zrecord(Iq:Iq)
 
C	    .Look for first letter of '.', 'F', or 'T'
          If (ThisChar .Ne. ' ') Then
 
             If (ThisChar .Eq. '.') Then
                 Iq = Iq+1
                 ThisChar = Zrecord(Iq:Iq) ! Take char after '.'
             End If   ! '.'
 
             If (ThisChar .Eq. 'T') Then
                 NLogical = .True.
             Else If (ThisChar .Eq. 'F') Then
                 NLogical = .False.
             Else
                 Ier = 1
      	   Go to 101
             End If
 
             Iqq =  Iq
C	       .Search until delimiter ' ' or ',' is found. Set pointer.
             Do 600 Ic=Iqq,Izlen
               Iq = Iq+1
      	 If (Zrecord(Iq:Iq) .Eq. ',' .Or.
     1               Zrecord(Iq:Iq) .Eq. ' '      ) Then
      	     Iptr = Iq+1
                   Go to 501
      		 End If    ! delimiter
  600        Continue
             Iptr = Izlen
 
             Go To 501
 
         End If  ! Not blank
 
  500     Continue
 
C	  .All blank Card -
        Ier = 1
        Go to 101
 
  501   Continue
 
      Return    ! End of Record
 
  101 Continue
      Return
      End
