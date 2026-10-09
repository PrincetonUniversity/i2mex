!--------------------------------------------------------------------
!
! Convert string to integer
!
FUNCTION atoi(string,i_in) RESULT (value)
  IMPLICIT NONE
  CHARACTER (len=*),intent(in) :: string
  INTEGER,intent(in) :: i_in
  INTEGER :: i, ii, max, sign, value
  CHARACTER (len=10), PARAMETER :: digit = '0123456789'
  i=i_in
  max = len(trim(string))
  value = 0
  sign = 1
  CALL skipbl(string,i)
  IF (string(i:i)=='-') THEN
    sign = -1
    i = i + 1
  END IF
  DO 
    ii = index(digit,string(i:i))
    if(ii==0) exit
    value = (value*10) + (ii-1)
    i=i+1
  END DO
  value = sign*value
END FUNCTION atoi

!--------------------------------------------------------------------
!
! Skip blanks in a string, used by ATOI
!
SUBROUTINE skipbl(string,i)
  IMPLICIT NONE
  CHARACTER(len=*) :: string
  INTEGER :: i, max

  max = len(string)

  IF (string(i:)==' ') THEN
    i = max
  ELSE
    DO i = i, max
      IF (string(i:i)/=' ') EXIT
    END DO
  END IF
  return
END SUBROUTINE skipbl
