subroutine tdb_nplas(d,inum)
  ! return number of non-impurity main plasma thermal species
  use trdatbuf_obj
  implicit NONE
  type (trdatbuf) :: d
  integer, intent(out) :: inum  ! number of species returned

  inum = d%ngmax_d

end subroutine tdb_nplas

subroutine tdb_nfast(d,inum)
  ! return number of fast ion species
  use trdatbuf_obj
  implicit NONE
  type (trdatbuf) :: d
  integer, intent(out) :: inum  ! number of species returned

  inum = d%nsfast_d

end subroutine tdb_nfast

subroutine tdb_plas(d,inum,zz,aa,ierr)
  ! return (in REAL*8 vectors) the Z & A values of main plasma species
  ! e.g. in a D,T plasma zz(1:2) = 1.0, 1.0 & aa(1:2) = 2.0, 3.0 returned.
  use trdatbuf_obj
  implicit NONE
  type (trdatbuf) :: d
  integer, intent(in) :: inum  ! number of species (see tdb_nplas)
  real*8, intent(out) :: zz(inum)  ! Z values (atomic numbers) returned.
  real*8, intent(out) :: aa(inum)  ! A values (atomic weights) returned.
  integer, intent(out) :: ierr ! error code (0=OK) (1 means inum too small)

  integer :: inumi,ia

  zz = 0
  aa = 0

  if(inum.lt.d%ngmax_d) then
     ierr = 1
  else
     inumi = d%ngmax_d
     ia = d%ladr_za
     zz(1:inumi) = d%datbuf(ia:ia+inumi-1)
     ia = ia + inumi
     aa(1:inumi) = d%datbuf(ia:ia+inumi-1)
  endif

end subroutine tdb_plas
      
subroutine tdb_fast(d,inum,zz,aa,ierr)
  ! return (in REAL*8 vectors) the Z & A values of fast ion species
  ! e.g. for D,He3 beams zz(1:2) = 1.0, 2.0 & aa(1:2) = 2.0, 3.0 returned.
  use trdatbuf_obj
  implicit NONE
  type (trdatbuf) :: d
  integer, intent(in) :: inum  ! number of species (see tdb_nfast)
  real*8, intent(out) :: zz(inum)  ! Z values (atomic numbers) returned.
  real*8, intent(out) :: aa(inum)  ! A values (atomic weights) returned.
  integer, intent(out) :: ierr ! error code (0=OK) (1 means inum too small)

  integer :: inumi,ia

  zz = 0
  aa = 0

  if(inum.lt.d%nsfast_d) then
     ierr = 1
  else
     inumi = d%nsfast_d
     ia = d%ladr_zaf
     zz(1:inumi) = d%datbuf(ia:ia+inumi-1)
     ia = ia + inumi
     aa(1:inumi) = d%datbuf(ia:ia+inumi-1)
  endif

end subroutine tdb_fast
