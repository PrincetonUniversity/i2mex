function r8_shafnl(Delprime)
  ! Returns the nonlinear Shafranov shift derivative "ShafNL",
  ! given the linear Shafranov shift derivative "Delprime".
  ! ShafNL is related to Delprime by the equation:
  !	ShafNL = Delprime * (1-ShafNL**2)**(3/2)
  ! This insures that -1 <= ShafNL <= 1.  Delprime may lie outside this
  ! interval, in which case the flux surfaces would nonsensically overlap
  ! if the nonlinear correction were not used.  SMP was used to solve the
  ! above cubic equation.
  use iso_c_binding, only: fp => c_double
  implicit none
  real(fp) :: r8_shafnl
  real(fp) :: delprime,b,xprs
  real(fp) :: third = (1._fp/3._fp)
 
  logical :: shafnl_flag
  data shafnl_flag /.true./ ! .false. = always use linear shift!
 
  ! For extremely small values of Delprime, the linear formula is perfectly
  ! O.K., while the nonlinear formula suffers from round-off errors.
 
  if(abs(Delprime) .le. 0.001_fp .or. .not. shafnl_flag) then
    r8_shafnl=Delprime
    return
  end if

  b=abs(Delprime)
  xprs=(-27.0_fp/b**2+(729.0_fp/b**4+108.0_fp/b**6)**0.5_fp)**third
  r8_shafnl=1.0_fp-2.0_fp**third/b**2/xprs + xprs/3/2**third

  r8_shafnl=sqrt(r8_shafnl)*SIGN (1.0_fp, Delprime)

  return
end function r8_shafnl 
