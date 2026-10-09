!******************** START FILE FPOLAR.FOR ; GROUP TKBNDRY ******************
!.....................................................................
!  REV DMC -- REMOVED TRANSP COMMON
!   SO I CAN USE THIS IN RPLOT
!
!  REV BB / MAY 94 -- USE ATAN2 INSTEAD OF ATAN TO EVALUATE FPOLAR
!
!   output in range [0,twopi] instead of range [-pi,pi]
!   atan2 and fpolar are the same in the upper half plane; fpolar
!   is twopi greater in the lower half plane.
!
function r8_fpolar ( ZR77, ZZ77 )
  !
  !	CALC. THE POLAR ANGLE DEFINED BY TAN(THETA)=ZZ77/ZR77
  !
  use iso_c_binding, only: fp => c_double
  use physconst_mod, only: zero, half, pi, twopi 
  implicit none
  real(fp) :: r8_fpolar
  real(fp) :: zz77, zr77
  !
  !--------------------------------------------
  !
  if (zz77 .ge. zero) then
    !
    !  upper half plane  --  including z=0 line
    !
    if (ZR77 .NE. zero) THEN
      !
      r8_fpolar = ATAN2 ( ZZ77 , ZR77 )
      !
    else
      !
      !  on the vertical axis, at or above the z=0 line
      !
      r8_fpolar = PI*half
      !
    end if
    !
  else
    !
    !  lower half plane  --  excluding z=0 line
    !
    r8_fpolar = ATAN2 ( ZZ77 , ZR77 ) + TWOPI
    !
  end if
  !
  return
end function r8_fpolar
