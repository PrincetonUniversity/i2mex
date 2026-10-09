module freeqbe_mod
 
  ! Free equilibrium format by Elena Belova
  ! pletzer@pppl.gov Thu Jun 28 16:03:33 EDT 2001
 
  integer, parameter :: feb_r8 = selected_real_kind(12,100)
  real(feb_r8), parameter :: feq_pi = 3.141592653589793116
 
  integer, parameter ::  feb_npsi0b=2, feb_nr0z0=2
 
  type freeqbe
 
     ! no of radial psi surfaces
     integer npsi
     ! grid resolution
     integer nr, nz
     ! units
     real(feb_r8) ulength_m, upsi_Wb__rad, up_Pa, uRB_mT
     ! magnetic axis position
     real(feb_r8) rzaxis(feb_nr0z0)
     ! psi axis-separatrix
     real(feb_r8) psi0b(feb_npsi0b)
 
     ! radial profiles
     real(feb_r8), dimension(:), pointer :: psi
     real(feb_r8), dimension(:), pointer :: p
     real(feb_r8), dimension(:), pointer :: RB
 
     ! grid/poloidal flux values
     real(feb_r8), dimension(:), pointer :: r, z
     real(feb_r8), dimension(:, :), pointer :: psi2
 
  end type freeqbe
 
contains
 
  subroutine feb_init(feb, filename, ier)
 
    use ezcdf
 
    implicit none
 
    type(freeqbe), intent(out) :: feb
    character*(*), intent(in) :: filename
    integer, intent(out) :: ier
 
    integer ncid, dimlens(3), iok
    character*4 xtype
 
    ier = 0
    call ezcdf_open(ncid, filename, 'r', iok)
    if (iok /= 0) then
       ier = 1
       return
    endif
 
    call cdfInqVar(ncid, 'psi', dimlens, xtype, iok)
    feb % npsi = dimlens(1)
    call cdfInqVar(ncid, 'r', dimlens, xtype, iok)
    feb % nr = dimlens(1)
    call cdfInqVar(ncid, 'z', dimlens, xtype, iok)
    feb % nz = dimlens(1)
 
    allocate(feb % psi(feb % npsi), stat=iok)
    allocate(feb % p(feb % npsi), stat=iok)
    allocate(feb % RB(feb % npsi), stat=iok)
    allocate(feb % r(feb % nr), stat=iok)
    allocate(feb % z(feb % nz), stat=iok)
    allocate(feb % psi2(feb % nr, feb % nz), stat=iok)
    if(iok /= 0) then
       ier = 2
    endif
 
    call cdfGetVar(ncid, 'ulength_m', feb % ulength_m, iok)
    call cdfGetVar(ncid, 'upsi_Wb__rad', feb % upsi_Wb__rad, iok)
    call cdfGetVar(ncid, 'up_Pa', feb % up_Pa, iok)
    call cdfGetVar(ncid, 'uRB_mT', feb % uRB_mT, iok)
    call cdfGetVar(ncid, 'rzaxis', feb % rzaxis, iok)
    call cdfGetVar(ncid, 'psi0b', feb % psi0b, iok)
    call cdfGetVar(ncid, 'psi', feb % psi, iok)
    call cdfGetVar(ncid, 'p', feb % p, iok)
    call cdfGetVar(ncid, 'RB', feb % RB, iok)
    call cdfGetVar(ncid, 'r', feb % r, iok)
    call cdfGetVar(ncid, 'z', feb % z, iok)
    call cdfGetVar(ncid, 'psi2', feb % psi2, iok)
    if(iok/=0) then
       ier = 3
    endif
 
    call ezcdf_close(ncid, iok)
 
  end subroutine feb_init
 
  subroutine feb_free(feb, ier)
    implicit none
    type(freeqbe), intent(inout) :: feb
    integer, intent(out) :: ier
 
    integer iok
 
    ier = 0
    deallocate(feb % psi, stat=iok)
    deallocate(feb % p, stat=iok)
    deallocate(feb % RB, stat=iok)
    deallocate(feb % r, stat=iok)
    deallocate(feb % z, stat=iok)
    deallocate(feb % psi2, stat=iok)
 
    feb % nr = 0
    feb % nz = 0
 
    if(iok /= 0) then
       ier = 4
    endif
  end subroutine feb_free
 
  subroutine feb_toMKS(feb, ier)
 
    ! Convert all data members to MKS units and
    ! set poloidal flux = 0 on axis.
 
    implicit none
    type(freeqbe), intent(inout) :: feb
    integer, intent(out) :: ier
 
    ier = 0
 
    feb% psi   = feb% psi  - feb% psi0b(1)
    feb% psi2  = feb% psi2 - feb% psi0b(1)
    feb% psi0b = feb% psi0b - feb% psi0b(1)
 
    feb% psi   = feb% psi  * feb% upsi_Wb__rad
    feb% psi2  = feb% psi2 * feb% upsi_Wb__rad
    feb% p     = feb% p    * feb% up_Pa
    feb% RB    = feb% RB   * feb% uRB_mT
    feb% r     = feb% r    * feb% ulength_m
    feb% z     = feb% z    * feb% ulength_m
 
    feb% rzaxis = feb% rzaxis * feb% ulength_m
    feb% psi0b  = feb% psi0b * feb% upsi_Wb__rad
 
    feb% upsi_Wb__rad = 1
    feb% up_Pa   = 1
    feb% uRB_mT   = 1
    feb% ulength_m   = 1
 
  end subroutine feb_toMKS
 
  subroutine feb_error(ier)
    implicit none
    integer, intent(in) :: ier
 
    select case (ier)
 
    case (1)
       print*,'freeqbe::freeqbe_init:: ERROR cannot open netCDF file!'
    case (2)
       print*,'freeqbe::freeqbe_init:: ERROR occurred while allocating!'
    case (3)
       print*,'freeqbe::freeqbe_init:: ERROR occurred while reading in!'
    case (4)
       print*,'freeqbe::freeqbe_free:: ERROR occurred while deallocating!'
 
    case default
 
    end select
  end subroutine feb_error
 
 
end module freeqbe_mod
 
