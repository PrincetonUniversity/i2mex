subroutine namels(iZ,iA,cname)
! names an element with charge and mass iZ and iA
!  (dmc 26 Apr 2001:  use Rob Andre's impurity labeler.  Note:  do not
!  real(fp) :: convert this routine!)

  use periodic_table_mod
  use iso_c_binding, only: fp => c_double

!============
!
  implicit none
  integer :: iZ  ! charge of the element
  integer :: iA  ! mass of the element
  real   r4A  ! atomic mass (real*4)  ! real*4 for periodic table software
  character cname*(*) ! character string name of the element

! default:  from periodic table...
  character(len=12) tmp
  
  integer :: iAmin,iZnum
  REAL   Astd  ! real*4 for periodic table software

! check: standard mass

  if(iZ .gt. 2) then

     iZnum=iZ
     do
        call standard_amu(iZnum,Astd)
        iAmin=Astd+1.0_fp
        if(iAmin.ge.iA) exit
        iZnum=iZnum+1
     end do

! impurity label
     
     r4A = iA
     tmp = to_periodic_table(iZnum,r4A,iZ,0)
     cname = tmp
     
  end if
  
! identify hydrogen isotopes:

  if(iZ .eq. 1)then
     if(iA .eq. 1) cname = 'H'
     if(iA .eq. 2) cname = 'D'
     if(iA .eq. 3) cname = 'T'
  end if
  
  if(iZ .eq. 2)then
     if(iA .eq. 3) cname = 'He3'
     if(iA .eq. 4) cname = 'He4'
  end if
  
  return
end subroutine namels
