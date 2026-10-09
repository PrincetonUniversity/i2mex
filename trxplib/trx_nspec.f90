subroutine trx_nspec(n_species)
!
!  return the number of plasma species (including electrons) present
!  in the currently open run.
!
!  if no run is open, return n_species=0.
!
  use trx_module
!
  implicit NONE
!
  integer, intent(out) :: n_species
!
!------------------------------------------------
!
  integer n_thi,n_thx,n_bi,n_rfi,n_fusi
  integer lunzer
!
!------------------------------------------------
!
  if(nsurf.eq.0) then
     write(lunzer(0),*) '?trx_nspec:  use trx_connect to open a run, first!'
     n_species=0
  else
     call rd_nspecies(n_species,n_thi,n_thx,n_bi,n_rfi,n_fusi)
  endif
end subroutine trx_nspec
