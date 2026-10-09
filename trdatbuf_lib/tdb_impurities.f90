!  The trdatbuf data may contain information of varying levels of detail on
!  plasma impurity content.
!
!  The traditional, simple model is to use a single "model impurity" known
!  as "tokamakium".  The charge number and atomic weight of this impurity
!  were allowed to be non-integers and vary in time, since "tokamakium"
!  represented a composite of some more complex mix of actual impurities
!  in a real tokamak.  If time variation is present, Zimp is a function of
!  time and the relation Aimp(t) = 2*Zimp(t) is assumed.
!
!  The more detailed model allows explicit specification of multi-species
!  impurity data.  In this case, there are multiple impurity elements and
!  also multiple charge states for each element, although the relative
!  concentrations of the various charge states may be estimated by a
!  "coronal equilibrium" approximation...

logical function tdb_tokamakium(d)

  ! return TRUE if the "tokamakium" model of impurity data is in effect

  use trdatbuf_obj
  type (trdatbuf) :: d

  tdb_tokamakium = ( d%NMIMP .eq. 0 )

end function tdb_tokamakium

!  multi impurity model:
!  return Z & A of impurity *elements*; see also tdb_imp_ions
!  vectors filled with zeroes if "tokamakium" model is in effect.

subroutine tdb_imp_elems(d,xzimps,aimps)
  use trdatbuf_obj
  use tdbsub_uts
  use periodic_table_mod
  implicit NONE
  INTEGER, PARAMETER  :: R8=SELECTED_REAL_KIND(12,100)
  type (trdatbuf) :: d
  real*8, dimension(:) :: xzimps,aimps  ! Z & A vectors
  !  fill in Z & A of impurity elements

  integer :: i,isize,imimps,lunmsg_tdb
  integer :: iz, ia,  ier
  real*8 :: r_mass, q_atom, q_ch

  if(d%nmimp.eq.0) return  ! tokamakium

  imimps = tdb_get_mimps(0)

  xzimps = ZERO
  aimps = ZERO

  isize=min(size(xzimps),size(aimps))
  if(isize.lt.imimps) then
     write(lunmsg_tdb(0),*) ' ?tdb_imp_elems: A & Z vectors too small:'
     write(lunmsg_tdb(0),*) '  size = ',isize,' expected = ',imimps
     return
  endif

  DO I = 1, imimps
     AIMPS(I)  = d%DATBUF(d%LAIMPS+I-1) ! load data -- assume good
     XZIMPS(I) = d%DATBUF(d%LXZIMPS+I-1)
     iz=XZIMPS(I) + 0.1_r8
     ia=AIMPS(I) + 0.5_r8
     call tr_species_convert(iz,iz,AIMPS(I),q_atom,q_ch,r_mass,ier)
     AIMPS(I)=r_mass/1.6726E-27_R8

  END DO

end subroutine tdb_imp_elems

!  multi impurity model:
!  return charge & atomic no. Z & and atomic weight A of impurity *ionss*; 
!  see also tdb_imp_elems
!  vectors filled with zeroes if "tokamakium" model is in effect.

subroutine tdb_imp_ions(d,xzimpx,xzimpxs,aimpx)
  use trdatbuf_obj
  use tdbsub_uts
  use periodic_table_mod
  implicit NONE
  INTEGER, PARAMETER  :: R8=SELECTED_REAL_KIND(12,100)
  type (trdatbuf) :: d
  real*8, dimension(:) :: xzimpx !  charge on ion e.g. Z(C+4)=4
  real*8, dimension(:) :: xzimpxs !  atomic number e.g. Z(C+4)=6
  real*8, dimension(:) :: aimpx  !  atomic weight e.g. A(C+4)=14 if Carbon-14

  !  fill in Z & A of impurity ions

  integer :: i,isize,imimpt,lunmsg_tdb
  integer :: iz, ia,  ier, iz_ch
  real*8 :: r_mass, q_atom, q_ch

  if(d%nmimp.eq.0) return  ! tokamakium

  imimpt = tdb_get_mimpt(0)

  xzimpx = ZERO
  xzimpxs = ZERO
  aimpx = ZERO

  isize=min(size(xzimpx),size(xzimpxs),size(aimpx))
  if(isize.lt.imimpt) then
     write(lunmsg_tdb(0),*) ' ?tdb_imp_ions: A & Z vectors too small:'
     write(lunmsg_tdb(0),*) '  size = ',isize,' expected = ',imimpt
     return
  endif
 
  DO I = 1, iMIMPT
     XZIMPX(I)  = d%DATBUF(d%LXZIMPX +I-1)
     XZIMPXS(I) = d%DATBUF(d%LXZIMPXS+I-1)
     AIMPX(I)   = d%DATBUF(d%LAIMPX  +I-1)
     iz=XZIMPXS(I) + 0.1_r8 !atomic #
     iz_ch=XZIMPXS(I) + 0.1_r8 !atomic charge of impurity
     ia=AIMPX(I) + 0.5_r8
     call tr_species_convert(iz,iz_ch,AIMPX(I),q_atom,q_ch,r_mass,ier)
     AIMPX(I)=r_mass/1.6726E-27_R8

  enddo

end subroutine tdb_imp_ions
