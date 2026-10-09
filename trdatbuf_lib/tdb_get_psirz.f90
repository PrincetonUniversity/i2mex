subroutine tdb_get_psirz(d,ztime,inumR,inumZ,PsiRZ,ier)

  !  get R and Z grids -- but, passed argument grid sizes must match or
  !  an error code is set and a message is printed.

  !  see sister subroutine "tdb_get_rzsizes"...

  use trdatbuf_obj
  use tdbsub_uts
  implicit NONE

  type (trdatbuf) :: d

  real*8, intent(in) :: ztime         ! time to which to interpolate
  integer,intent(in) :: inumR,inumZ   ! the grid sizes
  real*8 :: PsiRZ(inumR,inumZ)        ! Psi(R,Z) returned if sizes are correct
  integer, intent(out) :: ier         ! status code: 0=OK

  !------------------------------------
  integer :: nonlin,lunmsg_tdb,iloct,inumt,it1,ir,iz,ia1,ia2
  integer :: iloc_psi
  real*8 :: zfrac1,zfrac2
  !------------------------------------

  ier = 0

  nonlin = lunmsg_tdb(0)

  if(inumR.ne.d%nRpsi) then
     ier = ier + 1
     write(nonlin,*) ' ?tdb_get_Psirz: R grid size mismatch.'
     write(nonlin,*) '  Correct size is: ',d%nRpsi,'; passed size is: ',inumR
  endif

  if(inumZ.ne.d%nZpsi) then
     ier = ier + 1
     write(nonlin,*) ' ?tdb_get_Psirz: Z grid size mismatch.'
     write(nonlin,*) '  Correct size is: ',d%nZpsi,'; passed size is: ',inumZ
  endif

  if(ier.ne.0) return

  if(min(inumR,inumZ).le.0) then
     write(nonlin,*) ' %tdb_get_psirz: [R,Z] grids not found.'
     return
  endif

  iloct = d%ltpsi
  inumt = d%ntpsi

  if((iloct.eq.0).or.(inumt.eq.0)) then
     write(nonlin,*) '  Psi(t,R,Z) timebase address: ',iloct
     write(nonlin,*) '             timebase size:    ',inumt
     write(nonlin,*) '  --> zero not expected.'
     ier=2
     return
  endif

  call tdbsub_lookup(d%datbuf(iloct:iloct+inumt-1),inumt,ztime,it1,zfrac2)

  zfrac1 = ONE-zfrac2

  iloc_psi = d%lfpsi

  do iz=1,inumZ
     do ir=1,inumR
        ia1 = iloc_psi + (iz-1)*inumR*inumt + (ir-1)*inumt + (it1-1)
        ia2 = ia1 + 1
        PsiRZ(ir,iz)=zfrac1*d%datbuf(ia1) + zfrac2*d%datbuf(ia2)
     enddo
  enddo

end subroutine tdb_get_psirz
