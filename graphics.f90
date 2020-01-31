  subroutine i2mex_2plotmtv(nt1, ns, the, psi, filename, ier)

    ! Dump i2mex_o data into the ascii file named filename for 
    ! postprocessing using the free plotmtv program. 
    ! ier error flag (ok if =0)

    use i2mex_mod
    implicit none
    integer, intent(in) :: nt1, ns
    real(i2mex_r8), intent(in) :: the(nt1), psi(ns)
    character*(*), intent(in) :: filename
    integer, intent(out) :: ier

    integer, parameter ::  iu=1
    integer, parameter ::  markertype=13, markercolor=4, linetype=1, linecolor=1 
    integer iok, i, j
    real(i2mex_r8) :: z1(ns), z2(nt1,ns), x(nt1,ns), z(nt1, ns)
    

    ier=0
    open(iu, file=filename, form='formatted')

    ! p

    call i2mex_getP(ns, psi, z1, iok)
    call i2mex_error(iok)

    write(iu,*)'$DATA=COLUMN'
    write(iu,*)'%toplabel = p'
    write(iu,*)'%xlabel = psi'
    write(iu,*)'%linetype = ',linetype
    write(iu,*)'%linecolor= ',linecolor
    write(iu,*)'%markertype=',markertype
    write(iu,*)'%markercolor=',markercolor
    write(iu,*)'x y'
    do i=1, ns
       write(iu,'(2e16.8)')  psi(i), z1(i)
    enddo
    write(iu,*)' '
    
    ! p-prime

    call i2mex_getPP(ns, psi, z1, iok)
    call i2mex_error(iok)

    write(iu,*)'$DATA=COLUMN'
    write(iu,*)'%toplabel = "d p /d psi"'
    write(iu,*)'%xlabel = psi'
    write(iu,*)'%linetype = ',linetype
    write(iu,*)'%linecolor= ',linecolor
    write(iu,*)'%markertype=',markertype
    write(iu,*)'%markercolor=',markercolor
    write(iu,*)'x y'
    do i=1, ns
       write(iu,'(2e16.8)')  psi(i), z1(i)
    enddo
    write(iu,*)' '
    
    ! q

    call i2mex_getQ(ns, psi, z1, iok)
    call i2mex_error(iok)

    write(iu,*)'$DATA=COLUMN'
    write(iu,*)'%toplabel = q'
    write(iu,*)'%xlabel = psi'
    write(iu,*)'%linetype = ',linetype
    write(iu,*)'%linecolor= ',linecolor
    write(iu,*)'%markertype=',markertype
    write(iu,*)'%markercolor=',markercolor
    write(iu,*)'x y'
    do i=1, ns
       write(iu,'(2e16.8)')  psi(i), z1(i)
    enddo
    write(iu,*)' '
    
    ! q-prime

    call i2mex_getQP(ns, psi, z1, iok)
    call i2mex_error(iok)

    write(iu,*)'$DATA=COLUMN'
    write(iu,*)'%toplabel = "d q /d psi"'
    write(iu,*)'%xlabel = psi'
    write(iu,*)'%linetype = ',linetype
    write(iu,*)'%linecolor= ',linecolor
    write(iu,*)'%markertype=',markertype
    write(iu,*)'%markercolor=',markercolor
    write(iu,*)'x y'
    do i=1, ns
       write(iu,'(2e16.8)')  psi(i), z1(i)
    enddo
    write(iu,*)' '
    
    ! g

    call i2mex_getG(ns, psi, z1, iok)
    call i2mex_error(iok)

    write(iu,*)'$DATA=COLUMN'
    write(iu,*)'%toplabel = g'
    write(iu,*)'%xlabel = psi'
    write(iu,*)'%linetype = ',linetype
    write(iu,*)'%linecolor= ',linecolor
    write(iu,*)'%markertype=',markertype
    write(iu,*)'%markercolor=',markercolor
    write(iu,*)'x y'
    do i=1, ns
       write(iu,'(2e16.8)')  psi(i), z1(i)
    enddo
    write(iu,*)' '
    
    ! g-prime

    call i2mex_getGP(ns, psi, z1, iok)
    call i2mex_error(iok)

    write(iu,*)'$DATA=COLUMN'
    write(iu,*)'%toplabel = "d g /d psi"'
    write(iu,*)'%xlabel = psi'
    write(iu,*)'%linetype = ',linetype
    write(iu,*)'%linecolor= ',linecolor
    write(iu,*)'%markertype=',markertype
    write(iu,*)'%markercolor=',markercolor
    write(iu,*)'x y'
    do i=1, ns
       write(iu,'(2e16.8)')  psi(i), z1(i)
    enddo
    write(iu,*)' '

    ! surface plots

    call i2mex_getX(nt1, ns, the, psi, x, iok)
    call i2mex_error(iok)
    call i2mex_getZ(nt1, ns, the, psi, z, iok)
    call i2mex_error(iok)
    
    
    
    ! mesh

    call i2mex_getGP(ns, psi, z1, iok)
    call i2mex_error(iok)

    write(iu,*)'$DATA=COLUMN'
    write(iu,*)'%toplabel = "Mesh"'
    write(iu,*)'%xlabel = X'
    write(iu,*)'%ylabel = Z'
    write(iu,*)'%linetype = ',linetype
    write(iu,*)'%linecolor= ',linecolor
    !write(iu,*)'%markertype=',markertype
    !write(iu,*)'%markercolor=',markercolor
    write(iu,*)'x y'
    do i=1, ns
       do j=1, nt1
          write(iu,'(2e16.8)')  x(j,i), z(j,i)
       enddo
    enddo
    write(iu,*)' '

!!$    ! Grad-Shafranov equation lhs
!!$
!!$    call i2mex_getDelStar(nt1, ns, the, psi, z2, iok)
!!$    call i2mex_error(iok)
!!$    write(iu,*)'$DATA=CONTCURVE'
!!$    write(iu,*)'%contfill '
!!$    write(iu,*)'%toplabel = "Laplacian-star psi"'
!!$    write(iu,*)'%xlabel = X'
!!$    write(iu,*)'%ylabel = Z'
!!$    !write(iu,'(3e16.8)')  x(1,1), z(1,1), z2(1,1) ! 1st surface can be degenerated
!!$    do i=2, ns
!!$       do j=1, nt1-1
!!$          write(iu,'(3e16.8)')  x(j,i), z(j,i), z2(j,i)
!!$       enddo
!!$    enddo
!!$    
!!$    
!!$    ! Grad-Shafranov equation rhs
!!$
!!$    call i2mex_getXJphi(nt1, ns, the, psi, z2, iok)
!!$    call i2mex_error(iok)
!!$    write(iu,*)'$DATA=CONTCURVE'
!!$    write(iu,*)'%contfill '
!!$    write(iu,*)'%toplabel = "R J-phi"'
!!$    write(iu,*)'%xlabel = X'
!!$    write(iu,*)'%ylabel = Z'
!!$    !write(iu,'(3e16.8)')  x(1,1), z(1,1), z2(1,1) ! 1st surface can be degenerated
!!$    do i=2, ns
!!$       do j=1, nt1-1
!!$          write(iu,'(3e16.8)')  x(j,i), z(j,i), z2(j,i)
!!$       enddo
!!$    enddo
!!$    
    
    ! Grad-Shafranov equation error

    call i2mex_getGsError(nt1, ns, the, psi, z2, iok)
    call i2mex_error(iok)
    write(iu,*)'$DATA=CONTCURVE'
    write(iu,*)'%contfill '
    write(iu,*)'%toplabel = "Grad-Shafranov equation error"'
    write(iu,*)'%xlabel = X'
    write(iu,*)'%ylabel = Z'
    !write(iu,'(3e16.8)')  x(1,1), z(1,1), z2(1,1) ! 1st surface can be degenerated
    do i=2, ns
       do j=1, nt1-1
          write(iu,'(3e16.8)')  x(j,i), z(j,i), z2(j,i)
       enddo
    enddo
    
    close(iu)


  end subroutine i2mex_2plotmtv
    
  subroutine i2mex_2Dx(nt1, ns, the, psi, array, filename, ier)

    ! Dump i2mex_o mesh and array data into the ascii file named filename for 
    ! postprocessing using Data Explorer

    use i2mex_mod
    implicit none
    integer, intent(in) :: nt1, ns
    real(i2mex_r8), intent(in) :: the(nt1), psi(ns)
    real(i2mex_r8), intent(in) :: array(nt1, ns)
    character*(*), intent(in) :: filename
    integer, intent(out) :: ier

    integer, parameter ::  iu=1
    integer i, j
    real(i2mex_r8), dimension(:,:), allocatable :: x, z
    
    ier=0

    ALLOCATE(x(nt1, ns), z(nt1, ns))
    call i2mex_getX(nt1, ns, the, psi, x, ier)
    call i2mex_error(ier)
    call i2mex_getZ(nt1, ns, the, psi, z, ier)
    call i2mex_error(ier)

    open(iu, file=filename, form='formatted')

    write(iu,*)'# Data Explorer file format'
    write(iu,*)'# the irregular positions '
    write(iu,*)'# x y '
    write(iu,'(a,i6,a)') &
         & 'object 1 class array type float rank 1 shape 2 items ', &
         & nt1*(ns-1),'  data follows'
    do j=2, ns
       do i=1,nt1
          write(iu,'(2e14.6)') x(i,j), z(i,j)
       enddo
    enddo

    DEALLOCATE(x, z)

    write(iu,*)'# The regular connections: '
    write(iu,*)'object 2 class gridconnections counts ', ns-1, ' ', nt1
    write(iu,*)'# The data on the grid nodes '
    write(iu,*)'object 3 class array type float rank 0 items ', nt1*(ns-1),' data follows'
    do j=2, ns
       do i=1,nt1
          write(iu,'(2e14.6)') array(i,j)
       enddo
    enddo
    write(iu,*)'attribute "dep" string "positions" '
    write(iu,*)'object "irregular positions regular connections" class field '
    write(iu,*)'component "positions" value 1 '
    write(iu,*)'component "connections" value 2 '
    write(iu,*)'component "data" value 3 '
    write(iu,*)'end '

    close(iu)

  end subroutine i2mex_2Dx

