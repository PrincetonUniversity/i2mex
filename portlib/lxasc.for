C-----------------------------------------------------------------------
C  LXASC -- STRING COMPARE UTILITY
C
C    THIS ROUTINE WAS CREATED AT THE REQUEST OF J.P. JERAL AT JET, TO
C  AID IN THE PORTABILITY OF TRANSP AND RELATED CODE WHICH INVOLVES
C  STRING COMPARISONS FOR LEXICAL ORDERING.
C
C    THE PROBLEM IS THAT THE IBM EBCDIC CHARACTER CODE IS NOT ORDERED,
C  AND THE USE OF OPERATORS ".GT.", etc., IN FORTRAN DO NOT HAVE THE
C  EXPECTED RESULTS ON THE IBM.
C
C    THIS ROUTINE DOES NOT GIVE A PORTABLE SOLUTION, BUT, IT SHOULD
C  BE "PLUG COMPATIBLE" WITH THE IBM VERSION OF LXASC DEVELOPED BY
C  JEAN-PAUL.
C
      logical function LXASC(str1,oper,str2)
C
      implicit NONE
      character*(*) str1,str2
      character*2 oper
C
C-----------------------------------------------
C
      IF(OPER.EQ.'GT') THEN
        LXASC = str1.GT.str2
      ELSE IF(OPER.EQ.'GE') THEN
        LXASC = str1.GE.str2
      ELSE IF(OPER.EQ.'EQ') THEN
        LXASC = str1.EQ.str2
      ELSE IF(OPER.EQ.'LE') THEN
        LXASC = str1.LE.str2
      ELSE IF(OPER.EQ.'LT') THEN
        LXASC = str1.LT.str2
      ELSE IF(OPER.EQ.'NE') THEN
        LXASC = str1.NE.str2
      ELSE
        LXASC = .FALSE.
        WRITE(6,'('' ?LXASC:  OPERATOR "'',A,''" NOT RECOGNIZED'')')
     >      OPER
      ENDIF
C
      RETURN
      END
