subroutine trx_spec_edge(maxn,n_species,sgas,srcy,ierr)

  ! retrieve the gas flow and recycling sources -- scalars for each thermal
  ! plasma species

  implicit NONE

  integer, intent(in) :: maxn           ! max no. of species expected
  integer, intent(out) :: n_species     ! actual no. of species found

  real*8, intent(out) :: sgas(maxn)     ! gasflow sources, N/sec
  real*8, intent(out) :: srcy(maxn)     ! recycling sources, N/sec

  !  these are effective sources, i.e. #ionizations/sec inside the
  !  plasma core (a smaller number than the species' neutral influxes)

  integer, intent(out) :: ierr          ! completion code, 0=OK

  !----------------------------------

  real :: zz(maxn),aa(maxn)
  integer :: izc(maxn)
  integer :: i,igot,iertmp

  character*10 agas(maxn),arcy(maxn)
  character*20 zuns_mks
  real*8 :: zval

  integer :: lunzer
  !----------------------------------

  ierr = 0
  sgas = 0
  srcy = 0

  call rd_th_scedg(maxn,agas,arcy,zz,aa,izc,igot)
  if(igot.eq.-1) then
     write(lunzer(0),*) ' ?trx_spec_edge: error in rd_th_scedg call.'
     ierr=1
     return
  endif

  do i=1,igot
     call trx_scal(agas(i),zuns_mks,zval,iertmp)
     if(iertmp.ne.0) then
        write(lunzer(0),*) ' ?trx_spec_edge: failed to acquire: ',agas(i)
        ierr=1
     endif
     sgas(i)=zval

     call trx_scal(arcy(i),zuns_mks,zval,iertmp)
     if(iertmp.ne.0) then
        write(lunzer(0),*) ' ?trx_spec_edge: failed to acquire: ',arcy(i)
        ierr=1
     endif
     srcy(i)=zval
  enddo

end subroutine trx_spec_edge
