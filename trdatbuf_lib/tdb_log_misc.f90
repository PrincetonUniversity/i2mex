logical function tdb_zeff_simple(d,ilev)
  ! return TRUE if Zeff can be had by a simple interpolation
  use trdatbuf_obj
  implicit NONE
  type (trdatbuf) :: d
  integer, intent(out) :: ilev   ! 1: Zeff=f(t); 2: Zeff=f(x,t); 3: other

  if((d%nmzeff.eq.1).or.(d%nmzeff.eq.2)) then
     tdb_zeff_simple=.TRUE.
     ilev = d%nmzeff
  else
     tdb_zeff_simple=.FALSE.
     ilev = 3
  endif

end function tdb_zeff_simple

logical function tdb_imp_single(d,xzimp,aimp)
  ! return TRUE if the TRANSP run has but a single impurity species
  use trdatbuf_obj
  implicit NONE
  type (trdatbuf) :: d
  real*8, intent(out) :: xzimp  ! Z of impurity (0 if time dependent)
  real*8, intent(out) :: aimp   ! A of impurity (0 if time dependent)

  if(d%nmimp.eq.0) then
     tdb_imp_single = .TRUE.
     xzimp = d%datbuf(8)     ! can be zero...
     aimp = d%datbuf(9)      ! can be zero...
  else
     tdb_imp_single = .FALSE.
     xzimp = 0
     aimp = 0
  endif

end function tdb_imp_single
