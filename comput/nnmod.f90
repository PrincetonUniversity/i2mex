! This is a Fortran 90 implementation of some of the core functionality
! of the Julia Neural Network package "FLUX".  It can evaluate dense layers
! and networks composed of such layers, but does not have the ability
! to calculate gradients or otherwise train them.
module ffluxnn
  use iso_c_binding, only: r32=>c_float, r64=>c_double
  implicit none

  integer, parameter :: ACT_NONE      =  0
  integer, parameter :: ACT_SIGMOID   =  1
  integer, parameter :: ACT_TANH      =  2
  integer, parameter :: ACT_FAST_TANH =  3
  integer, parameter :: ACT_ELU       =  4

  type denselayer32
     real(r32), dimension(:,:), allocatable :: W  ! Weight matrix
     real(r32), dimension(:), allocatable :: b    ! bias vector
     integer :: nin, nout                         ! dimensions
     integer :: act                               ! activation fn index
  end type denselayer32

  type denselayer64
     real(r64), dimension(:,:), allocatable :: W  ! Weight matrix
     real(r64), dimension(:), allocatable :: b    ! bias vector
     integer :: nin, nout                         ! dimensions
     integer :: act                               ! activation fn index
  end type denselayer64

  type model32
     type(denselayer32), dimension(:), allocatable :: layer
     real(r32), dimension(:,:), allocatable        :: work
     integer                                       :: nlayers, ndata, nlabels
  end type model32

  type model64
     type(denselayer64), dimension(:), allocatable :: layer
     real(r64), dimension(:,:), allocatable        :: work
     integer                                       :: nlayers, ndata, nlabels
  end type model64

  interface newlayer
     module procedure randlayer32, randlayer64, datalayer32, datalayer64
  end interface newlayer

  interface freelayer
     module procedure freedense32, freedense64
  end interface freelayer

  interface layereval
     module procedure layereval32, layereval64
  end interface layereval

  interface showlayer
     module procedure showdense32, showdense64
  end interface showlayer

  interface freemodel
     module procedure freemodel32, freemodel64
  end interface freemodel

  interface modeleval
     module procedure modeleval32, modeleval64
  end interface modeleval

contains

  subroutine randlayer32(layer, nin, nout, act)
    type(denselayer32), intent(inout) :: layer
    integer, intent(in)               :: nin, nout
    integer, intent(in), optional     :: act

    if (allocated(layer%W)) deallocate(layer%W)
    if (allocated(layer%b)) deallocate(layer%b)
    layer%nin = nin;  layer%nout = nout
    layer%act = 0

    if ((nin.lt.1).or.(nout.lt.1)) return

    allocate(layer%W(nout,nin), layer%b(nout))

    call RANDOM_NUMBER(layer%W)
    layer%W = 2.0_r32*layer%W - 1.0_r32
    layer%b = 0.
    if (present(act)) layer%act = act
  end subroutine randlayer32

  subroutine randlayer64(layer, nin, nout, act)
    type(denselayer64), intent(inout) :: layer
    integer, intent(in)               :: nin, nout
    integer, intent(in), optional     :: act

    if (allocated(layer%W)) deallocate(layer%W)
    if (allocated(layer%b)) deallocate(layer%b)
    layer%nin = nin;  layer%nout = nout
    layer%act = 0

    if ((nin.lt.1).or.(nout.lt.1)) return

    allocate(layer%W(nout,nin), layer%b(nout))

    call RANDOM_NUMBER(layer%W)
    layer%W = 2.0_r64*layer%W - 1.0_r64
    layer%b = 0.
    if (present(act)) layer%act = act
  end subroutine randlayer64

  subroutine datalayer32(layer, weights, bias, nin, nout, act)
    type(denselayer32), intent(inout)          :: layer
    integer, intent(in)                        :: nin, nout, act
    real(r32), dimension(nout,nin), intent(in) :: weights
    real(r32), dimension(nout), intent(in)     :: bias

    if (allocated(layer%W)) deallocate(layer%W)
    if (allocated(layer%b)) deallocate(layer%b)
    layer%nin = nin;  layer%nout = nout
    layer%act = 0

    if ((nin.lt.1).or.(nout.lt.1)) return

    allocate(layer%W(nout,nin), layer%b(nout))

    layer%W = weights
    layer%b = bias
    layer%act = act
  end subroutine datalayer32

  subroutine datalayer64(layer, weights, bias, nin, nout, act)
    type(denselayer64), intent(inout)          :: layer
    integer, intent(in)                        :: nin, nout, act
    real(r64), dimension(nout,nin), intent(in) :: weights
    real(r64), dimension(nout), intent(in)     :: bias

    if (allocated(layer%W)) deallocate(layer%W)
    if (allocated(layer%b)) deallocate(layer%b)
    layer%nin = nin;  layer%nout = nout
    layer%act = 0

    if ((nin.lt.1).or.(nout.lt.1)) return

    allocate(layer%W(nout,nin), layer%b(nout))

    layer%W = weights
    layer%b = bias
    layer%act = act
  end subroutine datalayer64

  subroutine freedense32(layer)
    type(denselayer32), intent(inout) :: layer

    if (allocated(layer%W)) deallocate(layer%W)
    if (allocated(layer%b)) deallocate(layer%b)
    layer%nin = 0;  layer%nout = 0;  layer%act = 0
  end subroutine freedense32

  subroutine freedense64(layer)
    type(denselayer64), intent(inout) :: layer

    if (allocated(layer%W)) deallocate(layer%W)
    if (allocated(layer%b)) deallocate(layer%b)
    layer%nin = 0;  layer%nout = 0;  layer%act = 0
  end subroutine freedense64

  subroutine layereval32(layer, invec, outvec)
    type(denselayer32), intent(in)       :: layer
    real(r32), dimension(*), intent(in)  :: invec
    real(r32), dimension(*), intent(out) :: outvec

    integer :: irow

    outvec(1:layer%nout) = MATMUL(layer%W, invec(1:layer%nin)) + layer%b

    if (layer%act.lt.1) return

    select case (layer%act)
    case (ACT_SIGMOID)
       do irow=1,layer%nout
          outvec(irow) = sigmoid32(outvec(irow))
       enddo
    case (ACT_TANH)
       do irow=1,layer%nout
          outvec(irow) = TANH(outvec(irow))
       enddo
    case (ACT_FAST_TANH)
       do irow=1,layer%nout
          outvec(irow) = fast_tanh32(outvec(irow))
       enddo       
    case (ACT_ELU)
       do irow=1,layer%nout
          outvec(irow) = elu32(outvec(irow))
       enddo
    end select
  end subroutine layereval32

  subroutine layereval64(layer, invec, outvec)
    type(denselayer64), intent(in)       :: layer
    real(r64), dimension(*), intent(in)  :: invec
    real(r64), dimension(*), intent(out) :: outvec

    integer :: irow

    outvec(1:layer%nout) = MATMUL(layer%W, invec(1:layer%nin)) + layer%b

    if (layer%act.lt.1) return

    select case (layer%act)
    case (ACT_SIGMOID)
       do irow=1,layer%nout
          outvec(irow) = sigmoid64(outvec(irow))
       enddo
    case (ACT_TANH)
       do irow=1,layer%nout
          outvec(irow) = TANH(outvec(irow))
       enddo
    case (ACT_FAST_TANH)
       do irow=1,layer%nout
          outvec(irow) = fast_tanh64(outvec(irow))
       enddo  
    case (ACT_ELU)
       do irow=1,layer%nout
          outvec(irow) = elu64(outvec(irow))
       enddo
    end select
  end subroutine layereval64

  subroutine showdense32(layer)
    type(denselayer32), intent(in) :: layer

    integer :: irow

    write(*,'(A,I7,A,I7,A)')'Dense 32-bit layer with',&
         layer%nin,' input(s),',layer%nout,' output(s)'

    if ((layer%nin.lt.1).or.(layer%nout.lt.1)) return

    write(*,*)
    write(*,'(A)') 'Weights:'
    do irow=1,layer%nout
       write(*,'(6f12.6)') layer%W(irow,:)
    enddo

    write(*,*)
    write(*,'(A)') 'Biases:'
    write(*,'(6f12.6)') layer%b

    write(*,*)
    write(*,'(A,I4)') 'Activation function ',layer%act
  end subroutine showdense32

  subroutine showdense64(layer)
    type(denselayer64), intent(in) :: layer

    integer :: irow

    write(*,'(A,I7,A,I7,A)')'Dense 64-bit layer with',&
         layer%nin,' input(s),',layer%nout,' output(s)'

    if ((layer%nin.lt.1).or.(layer%nout.lt.1)) return

    write(*,*)
    write(*,'(A)') 'Weights:'
    do irow=1,layer%nout
       write(*,'(3f18.12)') layer%W(irow,:)
    enddo

    write(*,*)
    write(*,'(A)') 'Biases:'
    write(*,'(3f18.12)') layer%b

    write(*,*)
    write(*,'(A,I4)') 'Activation function ',layer%act
  end subroutine showdense64

  subroutine freemodel32(model)
    type(model32), intent(inout) :: model

    integer :: ilayer

    if (allocated(model%layer)) then
       do ilayer=1,size(model%layer)
          call freelayer(model%layer(ilayer))
       enddo
       deallocate(model%layer)
    endif
    if (allocated(model%work)) deallocate(model%work)
    model%nlayers = 0;  model%ndata = 0;  model%nlabels = 0
  end subroutine freemodel32

  subroutine freemodel64(model)
    type(model64), intent(inout) :: model

    integer :: ilayer

    if (allocated(model%layer)) then
       do ilayer=1,size(model%layer)
          call freelayer(model%layer(ilayer))
       enddo
       deallocate(model%layer)
    endif
    if (allocated(model%work)) deallocate(model%work)
    model%nlayers = 0;  model%ndata = 0;  model%nlabels = 0
  end subroutine freemodel64

  subroutine modeleval32(model, data, label)
    type(model32), intent(inout)         :: model
    real(r32), dimension(*), intent(in)  :: data
    real(r32), dimension(*), intent(out) :: label

    integer :: ilayer, jw=1

    if (model%nlayers.eq.1) then
       call layereval(model%layer(1), data, label)
    else
       call layereval(model%layer(1), data, model%work(:,jw))
       do ilayer=2,model%nlayers-1
          call layereval(model%layer(ilayer), model%work(:,jw), &
               model%work(:,3-jw))
          jw = 3 - jw
       enddo
       call layereval(model%layer(model%nlayers), model%work(:,jw), label)
    endif
  end subroutine modeleval32

  subroutine modeleval64(model, data, label)
    type(model64), intent(inout)         :: model
    real(r64), dimension(*), intent(in)  :: data
    real(r64), dimension(*), intent(out) :: label

    integer :: ilayer, jw=1

    if (model%nlayers.eq.1) then
       call layereval(model%layer(1), data, label)
    else
       call layereval(model%layer(1), data, model%work(:,jw))
       do ilayer=2,model%nlayers-1
          call layereval(model%layer(ilayer), model%work(:,jw), &
               model%work(:,3-jw))
          jw = 3 - jw
       enddo
       call layereval(model%layer(model%nlayers), model%work(:,jw), label)
    endif
  end subroutine modeleval64

  real(r32) function sigmoid32(x)
    real(r32), intent(in) :: x
    real(r32) :: t
     t = exp(-abs(x))
     if (x.ge.0.) then
        sigmoid32 = 1.0_r32/(1.0_r32 + t)
     else
        sigmoid32 = t/(1.0_r32 + t)
     endif

     ! ds/dx = s*(1-s)
  end function sigmoid32

  real(r64) function sigmoid64(x)
    real(r64), intent(in) :: x
    real(r64) :: t
     t = exp(-abs(x))
     if (x.ge.0.) then
        sigmoid64 = 1.0_r64/(1.0_r64 + t)
     else
        sigmoid64 = t/(1.0_r64 + t)
     endif
  end function sigmoid64

  real(r32) function fast_tanh32(x)
    real(r32), intent(in) :: x
    real(r32) :: t
    t = exp(x+x)
    fast_tanh32 = (t - 1.0_r32)/(t + 1.0_r32)

    ! df/dx = 4*t/(t+1)**2
  end function fast_tanh32

  real(r64) function fast_tanh64(x)
    real(r64), intent(in) :: x
    real(r64) :: t
    t = exp(x+x)
    fast_tanh64 = (t - 1.0_r64)/(t + 1.0_r64)
  end function fast_tanh64

  real(r32) function elu32(x)
    real(r32), intent(in) :: x
    if (x.ge.0.) then
       elu32 = x
    else
       elu32 = exp(x) - 1.0_r32
    endif
  end function elu32

  real(r64) function elu64(x)
    real(r64), intent(in) :: x
    if (x.ge.0.) then
       elu64 = x
    else
       elu64 = exp(x) - 1.0_r64
    endif
  end function elu64

end module ffluxnn
