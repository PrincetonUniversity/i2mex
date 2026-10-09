!
! --- getLog2Int1 ---
!
! fetch the integer 2**k+1>=ix for the smallest integer k>=2
!
function getLog2Int1(ix) result (kk)
  integer,intent(in) :: ix    ! target value
  
  integer :: kk    ! 2**k+1
  
  kk = 5
  do while(kk<ix)
     kk = 2*kk-1
  end do
end function getLog2Int1
