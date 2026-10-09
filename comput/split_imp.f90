subroutine split_imp(zimpm1,zimpm2,izimp1,izimp2,f1,f2,ierr)

  !  split an "average Z" impurity (Zimp) into two integer-Z impurities
  !  (izimp1), (izimp2) such that 

  !     f1*izimp1 +    f2*izimp2    = <Z> = zimpm1
  !     f1*izimp1**2 + f2*izimp2**2 = <Z**2> = zimpm2

  !  if:  izimp1 < min(<Z>,<Z**2>) < max(<Z>,<Z**2>) < izimp2  it can be
  !  shown that the fractions will be positive

  !  ierr=1 is set if min(f1,f2) = 0
  !  ierr=2 is set if min(f1,f2) < 0
  !  ierr=3 is set if izimp1=izimp2
  !  ierr=4 if min(iZimp1,iZimp2) < 1

  !  answer is:
  !     f1 = [<Z**2>-izimp2*<Z>]/[iZimp1*(iZimp1-iZimp2)]
  !     f2 = [<Z**2>-izimp2*<Z>]/[iZimp2*(iZimp2-iZimp1)]

  !  This can be used to create a pair of densities n1,n2 each associated
  !  with integer Z values, such that an "amalgum" of impurities with total
  !  density nimp and satisfying sum(Zx*nx)/nimp = <Z> and
  !  sum(Zx**2*nx)/nimp = <Z**2> can be represented as two densities n1,n2
  !  with Z values iZimp1 & iZimp2 respectively, and, contribution to Zeff 
  !  and quasineutrality are conserved;  n1=f1*nimp; n2=f2*nimp.
  !
  !  In the arguments of this routine the impurity density (nimp) is 
  !  normalized out; f1 and f2 are dimensionless.

  use iso_c_binding, only: fp => c_double
  use physconst_mod, only: zero
  implicit none

  !--------------------------------------
  real(fp), intent(in) :: Zimpm1,Zimpm2
  integer, intent(in) :: iZimp1,iZimp2

  real(fp), intent(out) :: f1,f2
  integer, intent(out) :: ierr
  !--------------------------------------

  ierr = 0

  f1 = ZERO
  f2 = ZERO

  if(izimp1.eq.izimp2) then
    ierr=3
    return
  endif

  if(min(iZimp1,iZimp2).lt.1) then
    ierr=4
    return
  endif

  f1 = (Zimpm2-iZimp2*Zimpm1)/(iZimp1*(iZimp1-iZimp2))
  f2 = (Zimpm2-iZimp1*Zimpm1)/(iZimp2*(iZimp2-iZimp1))

  if(min(f1,f2).le.ZERO) then
    if(min(f1,f2).eq.ZERO) then
      ierr=1
    else
      ierr=2
    endif
  endif

end subroutine split_imp
