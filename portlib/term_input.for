      SUBROUTINE TERM_input(jchar)
      use iso_c_binding, only: c_char
      implicit none
C
C
C...This routine was written by Harry H. Towner of the Princeton
C...Plasma Physics Lab.
C
C   Get a character from the terminal
C   --> DMC Oct 1991 <-- UNIX code using ordinary FORTRAN i/o
C                        which requires a carriage return
C
C...Parameters:
C*..jchar	- The byte that was input.
C
C***********************************************************************
C
      character(kind=c_char) jchar
C
#if __UNIX
      character(kind=c_char) :: ibuf
      character(kind=c_char) :: zgetc
 
      ibuf = zgetc()
      jchar=ibuf
#endif
C
      RETURN
      END
C--------------------------------
C
C  read in a single character from the terminal
C
      subroutine term_char_in(achar)
      use iso_c_binding, only: c_char
      implicit none
      character*1 achar
C
      character(kind=c_char) jchar
C
      call term_input(jchar)
      achar=jchar
C
      return
      end
 
