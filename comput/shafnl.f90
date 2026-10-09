      real function shafnl(Delprime)
! Returns the nonlinear Shafranov shift derivative "ShafNL",
! given the linear Shafranov shift derivative "Delprime".
! ShafNL is related to Delprime by the equation:
!	ShafNL = Delprime * (1-ShafNL**2)**(3/2)
! This insures that -1 <= ShafNL <= 1.  Delprime may lie outside this
! interval, in which case the flux surfaces would nonsensically overlap
! if the nonlinear correction were not used.  SMP was used to solve the
! above cubic equation.

      logical :: shafnl_flag
      data shafnl_flag /.true./ ! .false. = always use linear shift!

! For extremely small values of Delprime, the linear formula is perfectly
! O.K., while the nonlinear formula suffers from round-off errors.
      if(abs(Delprime) .le. .01 .or. .not. shafnl_flag) then
      shafnl=Delprime
      return
      end if

      b=abs(Delprime)
      xprs=(-27.0/b**2+(729.0/b**4+108.0/b**6)**0.5)**(1./3.)
      shafnl=1.0-2.0**(1./3.)/b**2/xprs + xprs/3/2**(1./3.)

      shafnl=Sqrt(shafnl)*SIGN (1.0, Delprime)
      return
      end


