C                                                 CHR_TOUPPER.FOR IN RPLOT_SUB
 
 
      CHARACTER FUNCTION CHR_TOUPPER (CHR)
C	                   CHR_TOUPPER CHANGES THE C*1 CHR TO UPPPER CASE IF
C                                      IT IS LOWER CASE.  TBT  8/90
      CHARACTER*1  CHR
 
C--------------
 
      IF (CHR .GE. 'a' .AND. CHR .LE. 'z')  THEN
 
          Jchr = Ichar('A') + (Ichar(Chr)-Ichar('a'))
          Chr  = Char(Jchr)
 
      END IF   ! CHR
 
      CHR_TOUPPER = CHR
 
      RETURN
      END
 
 
