integer function LKUPR_R8(X,TABLE,N)
  !
  ! Look up a real number X in a TABLE of length N.
  ! returns the value lkupr_r8 such that:
  !
  !	table(lkupr_r8) <= x < table(lkupr_r8+1)
  !
  ! assumes that table is in increasing order and does a binary search.
  !
  !============
  use iso_c_binding, only: fp => c_double
  implicit none
  integer, intent(in) :: n
  real(fp), intent(in) :: x
  real(fp), dimension(n), intent(in) :: table
  integer :: n1,n2,n3
  !
  if(X .LT. TABLE(1)) then
    LKUPR_R8=0
    GO TO 70000		!done
  end if
  !
  if(X .GE. TABLE(N)) then
    LKUPR_R8=N
    GO TO 70000		!done
  end if
  !
  N1=1
  N2=N
  !
100 continue	!LOOK AGAIN LOOP
  !
  N3=(N2+N1)/2
  !
  if(X .GE. TABLE(N3)) then
    N1=N3
  else
    N2=N3
  end if
  !
  if(N2 .EQ. N1+1) then
    LKUPR_R8=N1
    GO TO 70000	    !done
  else
    GOTO 100	    !LOOK AGAIN
  end if
 
70000 continue	!All returns from here
  !
  return
end function LKUPR_R8
