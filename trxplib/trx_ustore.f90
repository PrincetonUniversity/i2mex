subroutine trx_ustore(ierr_ck,ident,ixid,iordr,zunits,istat)
!
  use trx_module
  implicit NONE
!
!   add units label to profile label array
!
  integer, intent(inout) :: ierr_ck       ! error status code
  integer, intent(in) :: ident            ! xplasma id of profile
  integer, intent(in) :: ixid             ! x axis id (for 1d objects)
  integer, intent(in) :: iordr            ! interp. order (for 1d objects)
  !  ixid=iordr=0 expected for scalars
  character(*), intent(in) :: zunits      ! units label
  integer, intent(in) :: istat            ! units conversion flag
!
!   if ierr_ck=1, there was an error setting up the profile; just exit.
!   if ierr_ck=0, then ierr_ck will be set to 1 if  ident  is out of range.
!
!   ident -- xplasma index, being reused in trx_module to index into
!   arrays:  prof_units(...), mks_status(...)
!
!   zunits -- the units label
!
!   istat --   =0: MKS conversion was OK; =1: MKS conversion failure
!                                             TRANSP units used.
!-----------------------------------
  integer lunzer
!-----------------------------------
!
  if(ierr_ck.ne.0) return
!
  if((ident.le.0).or.(ident.gt.maxid)) then
     write(lunzer(0),*) ' ?trx_ustore:  item id = ',ident,' out of range.'
     ierr_ck=1
     return
  endif
!
  prof_units(ident)=zunits
  units_status(ident)=istat
  xaxis_id(ident)=ixid
  prof_iord(ident)=iordr
!
  call eqm_set_units(ident,zunits,ierr_ck)
!
  return
!
end subroutine trx_ustore
