!------------------------------------------------------------------
!  ST40YR  GENERATE SHOT YEAR FROM SHOT NUMBER BY HISTORICAL ANALYSIS
!
subroutine ST40YR(NSHOT,ZYEAR)
  implicit none
  !
  integer, intent(in) :: NSHOT
  character(len=*), intent(out) :: ZYEAR
  !
  integer :: ilen
  integer, dimension(:), allocatable :: idigs
  integer :: num_digits, ix, rem
  !
  num_digits = FLOOR(LOG10(REAL(NSHOT)+1))
  ALLOCATE(idigs(num_digits))
  rem = NSHOT
  do ix = 1, num_digits
    idigs(ix) = rem - (rem/10)*10 ! Take advantage of integer division
    rem = rem/10
  end do
  write(zyear,'(i2)') idigs(num_digits)*10+idigs(num_digits-1)
  deallocate(idigs)
  !
  return
end subroutine ST40YR
