!
! defines a single overloaded function for computing an integer hash from input data.  The
! floating point hashes include the first 12 digits or so and the exponent.
! 
!  hashCode(f,ihash)
!     f     - integer or real(fp) scalar or array data to hash
!     ihash - hash code input and output
!
module plasma_hash_code
  use iso_c_binding, only: fp => c_double
  use physconst_mod, only: zero, one
  implicit none

  interface hashCode
    module procedure hashCodeF0,hashCodeF1,hashCodeF2,hashCodeF3,hashCodeF4
    module procedure hashCodeI0,hashCodeI1,hashCodeI2,hashCodeI3,hashCodeI4
    module procedure hashCodeCn0,hashCodeCn1
    module procedure hashCdblF0,hashCdblF1,hashCdblF2,hashCdblF3,hashCdblF4
    module procedure hashCdblI0,hashCdblI1,hashCdblI2,hashCdblI3,hashCdblI4
    module procedure hashCdblCn0,hashCdblCn1
  end interface hashCode

  integer, parameter :: imax = 123456789  ! keep hash below this maximum
  private
  public hashCode

contains

  subroutine hashCodeCn0(cstr,k)
    character(len=*), intent(in) :: cstr
    integer,intent(inout) :: k   ! hash code

    integer :: ilen,i
    integer, dimension(:), allocatable :: iwk

    ilen = len(trim(cstr))
    allocate(iwk(ilen+1))

    iwk(1)=-ilen
    do i=1,ilen
      iwk(1+i)=ichar(cstr(i:i))
    enddo

    call hashCode_I32(ilen+1,iwk,k)
    deallocate(iwk)
  end subroutine hashCodeCn0

  subroutine hashCodeCn1(astr,k)
    character(len=*), intent(in), dimension(:) :: astr
    integer,intent(inout) :: k   ! hash code

    integer :: ilen,ilentot,isize,i,j
    integer, dimension(:), allocatable :: iwk

    ilentot = 0
    isize=size(astr)
    do j=1,isize
      ilen = len(trim(astr(j)))
      ilentot = ilentot + ilen + 1
    enddo

    allocate(iwk(ilentot))

    ilentot = 0
    do j=1,isize
      ilen = len(trim(astr(j)))
      iwk(ilentot+1) = -ilen
      do i=1,ilen
        iwk(ilentot+1+i)=ichar(astr(j)(i:i))
      enddo
    enddo

    call hashCode_I32(ilen+1,iwk,k)
    deallocate(iwk)
  end subroutine hashCodeCn1

  subroutine hashCodeF0(f,k)
    real(fp), intent(in)    :: f   ! data
    integer,intent(inout) :: k   ! hash code

    call hashCode_F64(1,(/f/),k)
  end subroutine hashCodeF0

  subroutine hashCodeF1(f,k)
    real(fp), intent(in)    :: f(:)   ! data
    integer,intent(inout) :: k      ! hash code

    call hashCode_F64(size(f),f,k)
  end subroutine hashCodeF1

  subroutine hashCodeF2(f,k)
    real(fp), intent(in)    :: f(:,:)   ! data
    integer,intent(inout) :: k        ! hash code
    integer :: n

    n = size(f,1)*size(f,2)

    call hashCode_F64(n,reshape(f,(/n/)),k)
  end subroutine hashCodeF2

  subroutine hashCodeF3(f,k)
    real(fp), intent(in)    :: f(:,:,:)   ! data
    integer,intent(inout) :: k          ! hash code
    integer :: n

    n = size(f,1)*size(f,2)*size(f,3)

    call hashCode_F64(n,reshape(f,(/n/)),k)
  end subroutine hashCodeF3

  subroutine hashCodeF4(f,k)
    real(fp), intent(in)    :: f(:,:,:,:)   ! data
    integer,intent(inout) :: k            ! hash code
    integer :: n

    n = size(f,1)*size(f,2)*size(f,3)*size(f,4)

    call hashCode_F64(n,reshape(f,(/n/)),k)
  end subroutine hashCodeF4

  subroutine hashCodeI0(f,k)
    integer,intent(in)    :: f   ! data
    integer,intent(inout) :: k   ! hash code

    call hashCode_I32(1,(/f/),k)
  end subroutine hashCodeI0

  subroutine hashCodeI1(f,k)
    integer,intent(in)    :: f(:)   ! data
    integer,intent(inout) :: k      ! hash code

    call hashCode_I32(size(f),f,k)
  end subroutine hashCodeI1

  subroutine hashCodeI2(f,k)
    integer,intent(in)    :: f(:,:)   ! data
    integer,intent(inout) :: k        ! hash code
    integer :: n

    n = size(f,1)*size(f,2)

    call hashCode_I32(n,reshape(f,(/n/)),k)
  end subroutine hashCodeI2

  subroutine hashCodeI3(f,k)
    integer,intent(in)    :: f(:,:,:)   ! data
    integer,intent(inout) :: k          ! hash code
    integer :: n

    n = size(f,1)*size(f,2)*size(f,3)

    call hashCode_I32(n,reshape(f,(/n/)),k)
  end subroutine hashCodeI3

  subroutine hashCodeI4(f,k)
    integer,intent(in)    :: f(:,:,:,:)   ! data
    integer,intent(inout) :: k            ! hash code
    integer :: n

    n = size(f,1)*size(f,2)*size(f,3)*size(f,4)

    call hashCode_I32(n,reshape(f,(/n/)),k)
  end subroutine hashCodeI4

  !-------------------------------------------

  subroutine hashCdblCn0(cstr,k1,k2)
    character(len=*), intent(in) :: cstr
    integer,intent(inout) :: k1,k2   ! hash code

    integer :: ilen,i
    integer, dimension(:), allocatable :: iwk

    ilen = len(trim(cstr))
    allocate(iwk(ilen+1))

    iwk(1)=-ilen
    do i=1,ilen
      iwk(1+i)=ichar(cstr(i:i))
    enddo

    call hashCdbl_I32(ilen+1,iwk,k1,k2)
    deallocate(iwk)
  end subroutine hashCdblCn0

  subroutine hashCdblCn1(astr,k1,k2)
    character(len=*), intent(in), dimension(:) :: astr
    integer,intent(inout) :: k1,k2   ! hash code

    integer :: ilen,ilentot,isize,i,j
    integer, dimension(:), allocatable :: iwk

    ilentot = 0
    isize=size(astr)
    do j=1,isize
      ilen = len(trim(astr(j)))
      ilentot = ilentot + ilen + 1
    enddo

    allocate(iwk(ilentot))

    ilentot = 0
    do j=1,isize
      ilen = len(trim(astr(j)))
      iwk(ilentot+1) = -ilen
      do i=1,ilen
        iwk(ilentot+1+i)=ichar(astr(j)(i:i))
      enddo
    enddo

    call hashCdbl_I32(ilen+1,iwk,k1,k2)
    deallocate(iwk)
  end subroutine hashCdblCn1

  subroutine hashCdblF0(f,k1,k2)
    real(fp), intent(in)    :: f      ! data
    integer,intent(inout) :: k1,k2  ! hash code

    call hashCdbl_F64(1,(/f/),k1,k2)
  end subroutine hashCdblF0

  subroutine hashCdblF1(f,k1,k2)
    real(fp), intent(in)    :: f(:)   ! data
    integer,intent(inout) :: k1,k2  ! hash code

    call hashCdbl_F64(size(f),f,k1,k2)
  end subroutine hashCdblF1

  subroutine hashCdblF2(f,k1,k2)
    real(fp), intent(in)    :: f(:,:)   ! data
    integer,intent(inout) :: k1,k2    ! hash code
    integer :: n

    n = size(f,1)*size(f,2)

    call hashCdbl_F64(n,reshape(f,(/n/)),k1,k2)
  end subroutine hashCdblF2

  subroutine hashCdblF3(f,k1,k2)
    real(fp), intent(in)    :: f(:,:,:)   ! data
    integer,intent(inout) :: k1,k2      ! hash code
    integer :: n

    n = size(f,1)*size(f,2)*size(f,3)

    call hashCdbl_F64(n,reshape(f,(/n/)),k1,k2)
  end subroutine hashCdblF3

  subroutine hashCdblF4(f,k1,k2)
    real(fp), intent(in)    :: f(:,:,:,:)   ! data
    integer,intent(inout) :: k1,k2        ! hash code
    integer :: n

    n = size(f,1)*size(f,2)*size(f,3)*size(f,4)

    call hashCdbl_F64(n,reshape(f,(/n/)),k1,k2)
  end subroutine hashCdblF4

  subroutine hashCdblI0(f,k1,k2)
    integer,intent(in)    :: f      ! data
    integer,intent(inout) :: k1,k2  ! hash code

    call hashCdbl_I32(1,(/f/),k1,k2)
  end subroutine hashCdblI0

  subroutine hashCdblI1(f,k1,k2)
    integer,intent(in)    :: f(:)   ! data
    integer,intent(inout) :: k1,k2  ! hash code

    call hashCdbl_I32(size(f),f,k1,k2)
  end subroutine hashCdblI1

  subroutine hashCdblI2(f,k1,k2)
    integer,intent(in)    :: f(:,:)   ! data
    integer,intent(inout) :: k1,k2    ! hash code
    integer :: n

    n = size(f,1)*size(f,2)

    call hashCdbl_I32(n,reshape(f,(/n/)),k1,k2)
  end subroutine hashCdblI2

  subroutine hashCdblI3(f,k1,k2)
    integer,intent(in)    :: f(:,:,:)   ! data
    integer,intent(inout) :: k1,k2      ! hash code
    integer :: n

    n = size(f,1)*size(f,2)*size(f,3)

    call hashCdbl_I32(n,reshape(f,(/n/)),k1,k2)
  end subroutine hashCdblI3

  subroutine hashCdblI4(f,k1,k2)
    integer,intent(in)    :: f(:,:,:,:)   ! data
    integer,intent(inout) :: k1,k2        ! hash code
    integer :: n

    n = size(f,1)*size(f,2)*size(f,3)*size(f,4)

    call hashCdbl_I32(n,reshape(f,(/n/)),k1,k2)
  end subroutine hashCdblI4

  !-------------------------------------------
  ! private routines...

  subroutine hashCode_F64(n, f, k)
    integer, intent(in)    :: n     ! number of elements
    real(fp),  intent(in)    :: f(:)  ! data
    integer, intent(inout) :: k     ! hash

    integer, parameter :: im8  = 100000000   

    integer :: ke, i8, i
    real(fp)  :: frac, f8

    k = mod(k,imax)
    do i=1, n
      ke = exponent(f(i))
      k  = mod(7*k+ke+imax,imax)

      frac = fraction(f(i))
      if (frac<zero) then
        k  = mod(7*k+1,imax)
        frac=abs(frac+one)
      end if

      f8 = frac*im8
      i8 = floor(f8)
      k  = mod(7*k+i8,imax)

      f8 = (f8-i8)*im8
      i8 = floor(f8)
      k  = mod(7*k+i8,imax)
      k  = abs(k)
    end do
  end subroutine hashCode_F64

  subroutine hashCode_I32(n, f, k)
    integer, intent(in)    :: n     ! number of elements
    integer, intent(in)    :: f(:)  ! data
    integer, intent(inout) :: k     ! hash

    integer, parameter :: im6  = 1000000   

    integer :: ke, i6, i
    integer :: frac, f6

    k = mod(k,imax)
    do i=1, n
      frac = f(i)
      if (frac<0) then
        k  = mod(7*k+1,imax)
        frac=abs(frac+imax)  ! make positive
      end if

      i6 = mod(frac,im6)      ! lower 6 digits
      k  = mod(7*k+i6,imax)

      i6 = frac/im6
      if (i6>0) then
        k = mod(7*k+i6,imax)
      end if
      k = abs(k)
    end do
  end subroutine hashCode_I32

  !-------------------------------------------

  subroutine hashCdbl_F64(n, f, k1,k2)
    integer, intent(in)    :: n     ! number of elements
    real(fp),  intent(in)    :: f(:)  ! data
    integer, intent(inout) :: k1,k2 ! hash

    integer, parameter :: im8  = 100000000   

    integer :: ke, i8, i, j
    real(fp)  :: frac, f8

    k1 = mod(k1,imax)
    k2 = mod(k2,imax)
    do i=1, n
      ke = exponent(f(i))
      k1  = mod(7*k1+ke+imax,imax)

      frac = fraction(f(i))
      if (frac<zero) then
        k1 = mod(7*k1+1,imax)
        frac=abs(frac+one)
      end if

      f8 = frac*im8
      i8 = floor(f8)
      k1 = mod(7*k1+i8,imax)

      f8 = (f8-i8)*im8
      i8 = floor(f8)
      k1 = mod(7*k1+i8,imax)
      k1 = abs(k1)

      j = n + 1 - i           ! scan in reverse order
      ke = exponent(f(j))
      k2  = mod(7*k2+ke+imax,imax)

      frac = fraction(f(j))
      if (frac<zero) then
        k2 = mod(7*k2+1,imax)
        frac=abs(frac+one)
      end if

      f8 = frac*im8
      i8 = floor(f8)
      k2 = mod(7*k2+i8,imax)

      f8 = (f8-i8)*im8
      i8 = floor(f8)
      k2 = mod(7*k2+i8,imax)
      k2 = abs(k2)
    end do
  end subroutine hashCdbl_F64

  subroutine hashCdbl_I32(n, f, k1,k2)
    integer, intent(in)    :: n     ! number of elements
    integer, intent(in)    :: f(:)  ! data
    integer, intent(inout) :: k1,k2 ! hash

    integer, parameter :: im6  = 1000000   

    integer :: ke, i6, i,j
    integer :: frac, f6

    k1 = mod(k1,imax)
    k2 = mod(k2,imax)
    do i=1, n
      frac = f(i)
      if (frac<0) then
        k1 = mod(7*k1+1,imax)
        frac=abs(frac+imax)  ! make positive
      end if

      i6 = mod(frac,im6)      ! lower 6 digits
      k1  = mod(7*k1+i6,imax)

      i6 = frac/im6
      if (i6>0) then
        k1 = mod(7*k1+i6,imax)
      end if
      k1 = abs(k1)

      j = n + 1 - i           ! scan in reverse order
      frac = f(j)
      if (frac<0) then
        k2 = mod(7*k2+1,imax)
        frac=abs(frac+imax)  ! make positive
      end if

      i6 = mod(frac,im6)      ! lower 6 digits
      k2  = mod(7*k2+i6,imax)

      i6 = frac/im6
      if (i6>0) then
        k2 = mod(7*k2+i6,imax)
      end if
      k2 = abs(k2)
    end do
  end subroutine hashCdbl_I32

end module plasma_hash_code
