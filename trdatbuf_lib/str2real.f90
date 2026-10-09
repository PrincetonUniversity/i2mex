!-----------------------------------------------------------------------
! Store an 8-byte string in an array of floats, or reverse the process
SUBROUTINE str2real(name, code, idir)
  USE ISO_C_BINDING, ONLY : r64 => C_DOUBLE
  IMPLICIT NONE
  INTRINSIC IACHAR, ACHAR

  ! Arguments
  character(len=8), intent(inout)        :: name  ! The string
  real(r64), dimension(8), intent(inout) :: code  ! The array
  integer, intent(in)                    :: idir  ! 1 => Encode;  -1 => Decode

  ! Local variable
  integer :: i

  if (idir.gt.0) then
     do i=1,8
        code(i) = real(iachar(name(i:i)),r64)
     enddo
  else
     do i=1,8
        name(i:i) = achar(int(code(i)))
     enddo
  endif

END SUBROUTINE str2real
