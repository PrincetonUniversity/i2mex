!
! ------------ get_fortran_const ----------
! Return the current fortran's concept of .true. and .false.
! For use by a C routine.
!
subroutine get_fortran_const(ftrue, ffalse)
  implicit none

  logical :: ftrue   ! fortran true
  logical :: ffalse  ! fortran false

  ftrue  = .true.
  ffalse = .false.
end subroutine get_fortran_const
