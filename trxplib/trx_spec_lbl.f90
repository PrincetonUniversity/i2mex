subroutine trx_spec_lbl(maxn,n_species,slbl,itype,Zr8,Ar8,IZc)
!
!  return the number of plasma species, and,
!  return a character*20 label for each plasma species present in the run, and,
!  return a (real*8) Z and A for each specie and IZc atomic number
!
!  for some model species, Z and A vary in time; therefore the time must
!  be selected before this routine is called.
!
  use trx_module
!
  implicit NONE
  integer, intent(in) :: maxn   ! max no. of species (dimension of slbl)
  integer, intent(out) :: n_species   ! actual no. of species found
  character*20, intent(out) :: slbl(maxn)   ! the species labels
  integer, intent(out) :: itype(maxn) ! type code of species (see below)
  real*8, intent(out) :: Zr8(maxn)    ! Z of species
  real*8, intent(out) :: Ar8(maxn)    ! A of species
  integer,intent(out) :: IZc(maxn)    ! atomic number of species
!
! return n_species=0 if an error occurs
!
! itype codes:  for j'th specie:
!
!   itype(j)=-1  -- electrons
!
!   itype(j)=+1  -- non-impurity thermal specie, usually H or He isotope
!                   can be Li
!   itype(j)=+2  -- impurity thermal specie:  Z and A are known, constant
!
!   itype(j)=+3  -- impurity thermal specie:  Z and A are functions of time
!                   "model impurity" -- could represent a hybrid; non-integer
!                   Z and A values possible.
!
!   itype(j)=+4  -- beam ion specie
!   itype(j)=+5  -- rf tail ion specie
!   itype(j)=+6  -- fusion product ion specie
!
!   NOTE ordering on output lists:
!     electrons
!     then non-impurity thermal ions
!     them impurity thermal ions
!     then non-thermal ions
!---------------------------------------
  character*10, dimension(:,:), allocatable :: abray
  integer, dimension(:), allocatable :: ifast
  real, dimension(:), allocatable :: aar4,zzr4
  integer i,iersum,ierr
  integer lunzer
  character*16 zunits
!---------------------------------------
!
  slbl=' '
  call trx_nspec(n_species)
  if(n_species.eq.0) return
  if(n_species.gt.maxn) then
     write(lunzer(0),*) '?trx_spec_lbl:  maxn=',maxn,' but: n_species=', &
          n_species
     n_species=0
     return
  endif
!
  if(itset.eq.0) then
     write(lunzer(0),*) '?trx_spec_lbl:  "time of interest" not yet selected.'
     n_species=0
     return
  endif
!
  allocate(abray(4,n_species))
  allocate(ifast(n_species))
  allocate(aar4(n_species),zzr4(n_species))
!
  call rd_species(n_species,n_species,slbl,abray,itype,ifast,zzr4,aar4,IZc)
  if(n_species.le.0) then
     n_species=0
  else
     iersum=0
     Zr8(1:n_species)=zzr4
     Ar8(1:n_species)=aar4
     do i=1,n_species
        if(itype(i).eq.3) then
           call trx_scal('XZIMP',zunits,Zr8(i),ierr)
           if(ierr.ne.0) iersum=iersum+1
           call trx_scal('AIMP',zunits,Ar8(i),ierr)
           if(ierr.ne.0) iersum=iersum+1
           IZc(i) = nint(Zr8(i))
        endif
     enddo
     if(iersum.gt.0) then
        write(lunzer(0),*) &
             '?trx_spec_lbl:  "AIMP" or "XZIMP" signals (f(t)) not found.'
        n_species=0
     endif
  endif
!
  deallocate(abray)
  deallocate(ifast)
  deallocate(aar4,zzr4)
!
end subroutine trx_spec_lbl
