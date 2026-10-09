      Subroutine TrCaps(Abbrev)
 
C	Change all lower case letters in Abbrev to Upper case.
 
 
      Character*(*) Abbrev
      Character*1   C
 
      Ilen = LEN(Abbrev)
 
      Do 100 I=1,Ilen
 
          C = Abbrev(I:I)
 
          If (C .Eq. 'a') Then
      	 C='A'
 
          Else if (C .Eq. 'b') Then
      	 C='B'
          Else if (C .Eq. 'c') Then
      	 C='C'
          Else if (C .Eq. 'd') Then
      	 C='D'
          Else if (C .Eq. 'e') Then
      	 C='E'
          Else if (C .Eq. 'f') Then
      	 C='F'
          Else if (C .Eq. 'g') Then
      	 C='G'
          Else if (C .Eq. 'h') Then
      	 C='H'
          Else if (C .Eq. 'i') Then
      	 C='I'
          Else if (C .Eq. 'j') Then
      	 C='J'
          Else if (C .Eq. 'k') Then
      	 C='K'
          Else if (C .Eq. 'l') Then
      	     C='L'
          Else if (C .Eq. 'm') Then
      	 C='M'
          Else if (C .Eq. 'n') Then
      	 C='N'
          Else if (C .Eq. 'o') Then
      	 C='O'
          Else if (C .Eq. 'p') Then
      	 C='P'
          Else if (C .Eq. 'q') Then
      	 C='Q'
          Else if (C .Eq. 'r') Then
      	 C='R'
          Else if (C .Eq. 's') Then
      	 C='S'
          Else if (C .Eq. 't') Then
      	 C='T'
          Else if (C .Eq. 'u') Then
      	 C='U'
          Else if (C .Eq. 'v') Then
      	 C='V'
          Else if (C .Eq. 'w') Then
      	 C='W'
          Else if (C .Eq. 'x') Then
      	 C='X'
          Else if (C .Eq. 'y') Then
      	 C='Y'
          Else if (C .Eq. 'z') Then
      	 C='Z'
          End If
 
          ABBREV(I:I) = C
 
  100   CONTINUE
 
      Return
      End
