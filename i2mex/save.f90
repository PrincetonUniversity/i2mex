  subroutine i2mex_save(filename, ier)
    
    ! Save i2mex_o object in netCDF file named filename.
    ! Use i2mex_load to load object from netCDF file.

    use i2mex_mod
    implicit none
    character*(*), intent(in) :: filename
    integer, intent(out) :: ier

    ier = 0
    call eq_save(filename,ier)
    if(ier/=0) ier = 50

  end subroutine i2mex_save

  subroutine i2mex_load(filename, ier)
    
    ! Load i2mex_o object from netCDF file named filename.

    use i2mex_mod
    implicit none
    character*(*), intent(in) :: filename
    integer, intent(out) :: ier

    integer iok, idum
    real(i2mex_r8), allocatable , dimension(:) :: psi

    ier = 0
    call eq_restore(filename,ier)
    if(ier/=0) ier = 51

    call eq_ganum('__RHO', i2mex_o%id_s)
    if(i2mex_o%id_s==0) then
       ier = 52
       print *,'failed to read id_s'
    endif
    call eq_ganum('__CHI', i2mex_o%id_t)
    if(i2mex_o%id_t==0) then
       ier = 53
       print *,'failed to read id_t'
    endif
    call eq_gfnum('PSI', i2mex_o%id_psi)
    if(i2mex_o%id_psi==0) then
       ier = 54
       print *,'failed to read id_psi'
    endif
    call eq_gfnum('P', i2mex_o%id_p)
    if(i2mex_o%id_p==0) then
       ier = 55
       print *,'failed to read id_p'
    endif
    call eq_gfnum('G', i2mex_o%id_g)
    if(i2mex_o%id_g==0) then
       ier = 56
       print *,'failed to read id_g'
    endif
    call eq_gfnum('Q', i2mex_o%id_q)
    if(i2mex_o%id_q==0) then
       ier = 57
       print *,'failed to read id_q'
    endif
    call eq_gfnum('R', i2mex_o%id_x)
    if(i2mex_o%id_x==0) then
       ier = 58
       print *,'failed to read id_x'
    endif
    call eq_gfnum('Z', i2mex_o%id_z)
    if(i2mex_o%id_z==0) then
       ier = 59
       print *,'failed to read id_z'
    endif

    call eq_ngrid(i2mex_o%id_s, i2mex_o%ns)
    if(i2mex_o%ns==0) then
       print *,'failed to get i2mex_o%ns'
    endif
    call eq_ngrid(i2mex_o%id_t, i2mex_o%nt1)
    if(i2mex_o%nt1==0) then
       print *,'failed to get i2mex_o%nt1'
    endif

!!$    allocate(psi(i2mex_o%ns))
!!$
!!$    call eq_grid(i2mex_o%id_psi,psi,i2mex_o%ns,idum,iok)
!!$    i2mex_o%psi_edge = psi(i2mex_o%ns)
!!$    i2mex_o%psi_axis = psi(1)
!!$
!!$    deallocate(psi)

  end subroutine i2mex_load

