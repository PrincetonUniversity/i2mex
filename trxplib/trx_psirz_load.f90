subroutine trx_psirz_load(ierr)

  use trx_module
  use xplasma_obj_instance
  implicit NONE

  ! load trxplib Psi(R,Z) (from TRANSP free bdy run) into xplasma2

  integer, intent(out) :: ierr

  !----------------------------------------------------
  ! local:

  integer :: id_psi,idprev,id_Rgrid,id_Zgrid,iertmp
  integer :: lunzer
  !----------------------------------------------------
  ierr = 0

  call xplasma_profID(s,'psi',id_psi)
  if(id_psi.eq.0) then
     ierr=1
     write(lunzer(0),*) ' ?trx_psirz_load: 1d psi(rho) not found.'
     return
  endif

  idprev = id_psi

  call xplasma_gridID(s,'__Rgrid',id_Rgrid)
  call xplasma_gridID(s,'__Zgrid',id_Zgrid)

  if(id_Rgrid.eq.0) then
     ierr=1
     write(lunzer(0),*) ' ?trx_psirz_load: __Rgrid ID: not found.'
  endif

  if(id_Zgrid.eq.0) then
     ierr=1
     write(lunzer(0),*) ' ?trx_psirz_load: __Zgrid ID: not found.'
  endif

  if(ierr.ne.0) return

  ! OK:

  call xplasma_author_set(s,'trxplib',iertmp)

  call xplasma_create_2dprof(s,'PSI_RZ', &
       id_Rgrid,id_Zgrid, PsiRZ_free, id_psi, ierr, &
       ispline=1, assoc_id=idprev, &
       label='Poloidal Flux',units='Wb/rad')

  call xplasma_author_clear(s,'trxplib',iertmp)

  if(ierr.ne.0) then
     write(lunzer(0),*) &
          ' ?trx_psirz_load: Psi(R,Z) xplasma_create_2dprof error:'
     call xplasma_error(s,ierr,lunzer(0))
  endif

end subroutine trx_psirz_load
