!
! C callable versions of the periodic table routines.  Only the pseudo-periodic table
! versions are implemented here.  The C++ declarations are in periodic_table.h
!
 
 
!
! ---------------------- cignore_user_amu ---------------
! set get periodic_table::ignore_user_AMU

subroutine cignore_user_amu(i)
  use periodic_table_mod
  implicit none

  integer :: i    ! -1 -> return ignore_user_AMU
                  !  0 -> set ignore_user_AMU=.false.
                  !  1 -> set ignore_user_AMU=.true.

  if (i>=0) then
     ignore_user_AMU = i>0
  end if

  if (ignore_user_AMU) then
     i=1
  else
     i=0
  end if
end subroutine cignore_user_amu

!
! ---------------------- cignore_electron_amu ---------------
! set get periodic_table::ignore_electron_AMU

subroutine cignore_electron_amu(i)
  use periodic_table_mod
  implicit none

  integer :: i    ! -1 -> return ignore_electron_AMU
                  !  0 -> set ignore_electron_AMU=.false.
                  !  1 -> set ignore_electron_AMU=.true.

  if (i>=0) then
     ignore_electron_AMU = i>0
  end if

  if (ignore_electron_AMU) then
     i=1
  else
     i=0
  end if
end subroutine cignore_electron_amu

!
! ---------------------- cps_params ---------------
! set get periodic_table::ps_params

subroutine cps_params(i)
  use periodic_table_mod
  implicit none

  integer :: i    ! -1 -> return ps_params
                  !  0 -> set ps_params=.false.
                  !  1 -> set ps_params=.true.

  if (i>=0) then
    ps_params = i>0
  end if

  if (ps_params) then
    i=1
  else
    i=0
  end if

  return
end subroutine cps_params

!
! --------------------- cto_periodic_table ---------------------
! return the periodic table symbol except the isotopes of hydrogen
! are returned as H,D,T if ipseudo/=0
!
! cout = null terminated char* with space for at least 12 characters+null
!
! extern "C" void cto_periodic_table(const int& ipseudo, const int& z, const float& a,
!                                    const int& charge,  const int& isymbol, char* cout)
!
subroutine cto_periodic_table(ipseudo, z, a, charge, isymbol, cout)
  use iso_c_binding, only: c_char
  use periodic_table_mod
  implicit none
 
  integer, intent(in) :: ipseudo ! 0 -> standard periodic table
                                 ! 1  -> return H,D,T for hydrogen isotopes
  integer, intent(in) :: z       ! the atomic number of the element
  real(fp), intent(in) :: a      ! the atomic weight of the element
  integer, intent(in) :: charge  ! the charge on the element
  integer, intent(in) :: isymbol ! symbol to use to separate element name and charge
                                 ! 0 -> _  this is an underscore
                                 ! 1 -> +
 
  character(kind=c_char), intent(out) :: cout(13) ! the null terminated symbol returned
                                        ! from to_pseudo_periodic_table
 
  character(len=12) :: cin              ! buffer for fortran symbol
 
  cin = to_pseudo_periodic_table(ipseudo, z, a, charge, isymbol)
 
  call cstring(cin,cout,'2C')

  return
end subroutine cto_periodic_table
 
 
!
! --------------------- cinv_periodic_table ---------------------
!
! return the inv_periodic_table data where the hydrogen isotopes
! are given by H,D,T if ipseudo/=0
!
! symbol = null terminated char* with no more then 12 characters
!
! extern "C" void cinv_periodic_table(const int& ipseudo, const char* symbol,
!                                     const int& idefault, int& z, float& a, int& charge)
!
subroutine cinv_periodic_table(ipseudo, symbol, idefault, z, a, charge)
  use iso_c_binding, only: c_char
  use periodic_table_mod
  implicit none
 
  integer, intent(in)      :: ipseudo   ! 0 -> standard periodic table
                                        ! 1 -> use H,D,T for hydrogen isotopes
  character(kind=c_char), intent(in) :: symbol    ! null terminated string containing element (e.g. 'C12+4 ')
                                        ! must be left justified
  integer, intent(in)      :: idefault  ! nonzero to return an atomic weight even
                                        ! if data is not contained in the string
 
  integer, intent(out) :: z             ! atomic number of the element
  real(fp), intent(out) :: a            ! atomic weight of the element
  integer, intent(out) :: charge        ! charge state of the element
 
  logical           :: default
  character(len=12) :: fsymbol
 
  default = (idefault/=0)
  call cstring(fsymbol,symbol,'2F')
  call inv_pseudo_periodic_table(ipseudo, fsymbol, default, z, a, charge)

  return
end subroutine cinv_periodic_table
 
!
! -------------------- cname_periodic_table --------------------
!
! return the name of the element with the atomic number if A>0,
! the hydrogen isotopes will be returned as their name if ipseudo>0
!
! cout = null terminated char* with space for at least 16 characters+null
!
! extern "C" void cname_periodic_table(const int& ipseudo, const int& z, const float& a, char* cout)
!
subroutine cname_periodic_table(ipseudo, z, a, cout)
  use iso_c_binding, only: c_char
  use periodic_table_mod
  implicit none
 
  integer, intent(in) :: ipseudo ! 0 -> standard periodic table
                                 ! 1  -> return H,D,T for hydrogen isotopes
  integer, intent(in) :: z       ! the atomic number of the element
  real(fp), intent(in) :: a      ! the atomic weight of the element
 
  character(kind=c_char), intent(out) :: cout(16) ! the null terminated symbol returned
                                        ! from name_pseudo_periodic_table
 
  character(len=16) :: cin       ! buffer for fortran name
 
  cin = name_pseudo_periodic_table(ipseudo, z, a)
 
  call cstring(cin,cout,'2C')

  return
end subroutine cname_periodic_table

!
! -------------------- cname_periodic_table --------------------
!
! return the name of the element with the atomic number if A>0,
! the hydrogen isotopes will be returned as their name if ipseudo>0
!
! cout = null terminated char* with space for at least 16 characters+null
!
! extern "C" void cname_periodic_table(const int& ipseudo, const int& z, const float& a, char* cout)
!
subroutine ctr_species_convert_r8(iZatom, iZcharge, zAMU, qatom, qcharge, mass, ierr) &
     bind(C, name="ctr_species_convert_r8")
  use periodic_table_mod
  use iso_c_binding, only:  c_int, fp => c_double
  implicit none
 
  integer :: iZatom    ! atomic number: 1 for H, 2 for He, ...

  integer :: iZcharge  ! ionic charge btw 1 and iZatom

  REAL(fp)  :: zAMU      ! atomic weight in AMU, nearest integer

  REAL(fp)  :: qatom     ! charge of fully stripped ion
  REAL(fp)  :: qcharge   ! actual charge
  REAL(fp)  :: mass      ! mass in kg

  integer :: ierr      ! status code, 0=OK

  integer :: i,j

  !
  !Used to create element_index() array
  !
  !do i=1, 109
  !   write(*,'(i3,a)',advance='no') j,', '
  !end do

  call tr_species_convert_r8(iZatom, iZcharge, zAMU, qatom, qcharge, mass,ierr)
end subroutine ctr_species_convert_r8
