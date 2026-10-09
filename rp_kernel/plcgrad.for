 
      Subroutine PLCGRAD( text, Nptr)
 
C     Called from PLCPARSE:
C     DIV operator found in text (Nptr points at next character after DIV)
 
C     Insert into text, such that
 
C     DIV(F) becomes
 
C     FLXDIFF(SURF*(F))/DVOL
 
      Implicit None
 
 
      Character*(*) text
      Character*170 TEMPTEXT
 
      Integer Nptr
 
      Integer Nend, I, J, Nparen
 
 
C     ------------------------------------------------
 
      Nend = LEN(text)
 
      Do 200 I=Nptr,Nend-1
 
         IF (text(I:I) .EQ. '(') THEN
 
            Nparen = 1
 
            Do 100 J=I+1,Nend
 
               IF (text(J:J) .EQ. ')' ) Nparen = Nparen-1
 
               IF (Nparen .EQ. 0) THEN
 
C                 .Found matching parenthesis - close GRAD
                  temptext = text(1:I)   // 'SURF*('      //
     1                       text(I+1:j) // ')/DVOL' //
     2                       text(J+1:Nend)
                  text     = temptext
                  text(Nend:Nend) = ';'     ! just in case other was wiped out.
                  GO TO 8000
               ENDIF
 
               IF (text(J:J) .EQ. '(' ) Nparen = Nparen+1     ! count balancing '(')'
 
 100        Continue  ! Do j=
 
            ! Error - should never reach here.
 
         ENDIF
 
 200  Continue     ! do i=
 
      ! Should never reach here!!
 
 8000 Continue
 
      Return
      End
