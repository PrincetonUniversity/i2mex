!
! submsg_mod
!
! Defines a submsg type which can be used for returning a multi-line
! error message from a Fortran subroutine.
!
! A default submsg is defined with C-callable accessor routines.
!

module submsg_mod

  implicit none

  integer, parameter :: nsubmsg = 132

  ! Single line for temporary use in external code.
  character(len=nsubmsg) :: mline

  !
  ! ------------------------------------------------------------------
  ! submsg type
  ! ------------------------------------------------------------------
  !

  type :: submsg
     private

     ! Number of logical lines.
     ! -1 means storage has not yet been allocated.
     integer :: num = -1

     ! Message storage.
     character(len=nsubmsg), pointer, dimension(:) :: msg => null()

  end type submsg


  !
  ! ------------------------------------------------------------------
  ! Default submsg
  ! ------------------------------------------------------------------
  !

  logical, private, save :: igot = .false.

  type(submsg), pointer, private, save :: default_msg => null()


  !
  ! ------------------------------------------------------------------
  ! Generic interfaces
  ! ------------------------------------------------------------------
  !

  interface add_msg
     module procedure add_msg_char
     module procedure add_msg_msg
  end interface add_msg


  interface copy_msg
     module procedure copy_msg_char
     module procedure copy_msg_msg
  end interface copy_msg


contains


  !
  ! ------------------------------------------------------------------
  ! new_msg
  !
  ! Initialize a new submsg.
  !
  ! The optional ns argument specifies an initial physical storage size.
  ! ------------------------------------------------------------------
  !

  subroutine new_msg(m, ns)

    type(submsg), intent(inout)   :: m
    integer, optional, intent(in) :: ns

    if (present(ns)) then
       if (ns > 0) then
          call resize_msg(m, ns)
       end if
    end if

  end subroutine new_msg


  !
  ! ------------------------------------------------------------------
  ! size_msg
  !
  ! Return the logical number of lines in a message.
  ! ------------------------------------------------------------------
  !

  function size_msg(m) result(n)

    type(submsg), intent(in) :: m
    integer                  :: n

    n = max(0, m%num)

  end function size_msg


  !
  ! ------------------------------------------------------------------
  ! resize_msg
  !
  ! ns >= 0:
  !     Set physical storage size to ns.
  !
  ! ns < 0:
  !     Ensure that at least abs(ns) free lines are available.
  ! ------------------------------------------------------------------
  !

  subroutine resize_msg(m, ns)

    type(submsg), intent(inout) :: m
    integer,      intent(in)    :: ns

    integer :: n
    integer :: k

    character(len=nsubmsg), pointer, dimension(:) :: p

    !
    ! Determine requested physical size.
    !

    if (ns >= 0) then

       n = ns

    else

       if (m%num <= 0) then

          n = abs(ns)

       else

          n = m%num + abs(ns)

          ! Already enough free space.
          if (associated(m%msg)) then
             if (n <= size(m%msg)) then
                return
             end if
          end if

          ! Increase with some extra space.
          n = max(n + 3, 4*n/3)

       end if

    end if


    !
    ! Allocate for the first time.
    !

    if (m%num < 0) then

       allocate(m%msg(max(0, n)))

       m%num = 0


    !
    ! Resize existing storage.
    !

    else

       if (associated(m%msg)) then
          if (n == size(m%msg)) then
             return
          end if
       end if

       allocate(p(max(0, n)))

       k = min(size(p), m%num)

       if (k > 0) then
          p(1:k) = m%msg(1:k)
       end if

       if (associated(m%msg)) then
          deallocate(m%msg)
       end if

       m%msg => p

       m%num = k

    end if

  end subroutine resize_msg


  !
  ! ------------------------------------------------------------------
  ! trim_msg
  !
  ! Remove ns logical lines from the end of the message.
  ! ------------------------------------------------------------------
  !

  subroutine trim_msg(m, ns)

    type(submsg), intent(inout) :: m
    integer,      intent(in)    :: ns

    if (m%num <= 0) then
       return
    end if

    if (ns <= 0) then
       return
    end if

    m%num = max(0, m%num - ns)

  end subroutine trim_msg


  !
  ! ------------------------------------------------------------------
  ! clear_msg
  !
  ! Empty the submsg and release storage.
  ! ------------------------------------------------------------------
  !

  subroutine clear_msg(m)

    type(submsg), intent(inout) :: m

    if (associated(m%msg)) then
       deallocate(m%msg)
       nullify(m%msg)
    end if

    m%num = -1

  end subroutine clear_msg


  !
  ! ------------------------------------------------------------------
  ! add_msg_char
  !
  ! Append one string to a message.
  ! ------------------------------------------------------------------
  !

  subroutine add_msg_char(m, c)

    type(submsg),      intent(inout) :: m
    character(len=*), intent(in)    :: c

    integer :: n

    call resize_msg(m, -1)

    n = m%num + 1

    m%msg(n) = c

    m%num = n

  end subroutine add_msg_char


  !
  ! ------------------------------------------------------------------
  ! add_msg_msg
  !
  ! Append all lines from one submsg to another.
  ! ------------------------------------------------------------------
  !

  subroutine add_msg_msg(m, mc)

    type(submsg), intent(inout) :: m
    type(submsg), intent(in)    :: mc

    integer :: i
    integer :: n

    if (mc%num <= 0) then
       return
    end if

    call resize_msg(m, -mc%num)

    do i = 1, mc%num

       n = m%num + 1

       m%msg(n) = mc%msg(i)

       m%num = n

    end do

  end subroutine add_msg_msg


  !
  ! ------------------------------------------------------------------
  ! copy_msg_char
  !
  ! Replace a message with one string.
  ! ------------------------------------------------------------------
  !

  subroutine copy_msg_char(m, c)

    type(submsg),      intent(inout) :: m
    character(len=*), intent(in)    :: c

    if (m%num <= 0) then
       call resize_msg(m, -1)
    end if

    m%msg(1) = c

    m%num = 1

  end subroutine copy_msg_char


  !
  ! ------------------------------------------------------------------
  ! copy_msg_msg
  !
  ! Replace one message with another message.
  ! ------------------------------------------------------------------
  !

  subroutine copy_msg_msg(m, mc)

    type(submsg), intent(inout) :: m
    type(submsg), intent(in)    :: mc

    integer :: i

    if (mc%num <= 0) then

       if (m%num > 0) then
          m%num = 0
       end if

       return

    end if


    if (m%num < mc%num) then
       call resize_msg(m, mc%num)
    end if


    do i = 1, mc%num
       m%msg(i) = mc%msg(i)
    end do

    m%num = mc%num

  end subroutine copy_msg_msg


  !
  ! ------------------------------------------------------------------
  ! line_msg
  !
  ! Return one line from a message.
  !
  ! Fortran indexing is used here: 1 ... size_msg(m).
  ! ------------------------------------------------------------------
  !

  subroutine line_msg(m, i, line)

    type(submsg),      intent(in)    :: m
    integer,           intent(in)    :: i
    character(len=*), intent(inout) :: line

    integer :: n

    character(len=nsubmsg) :: errmsg

    n = max(0, m%num)

    if ((i <= 0) .or. (i > n)) then

       write(errmsg, '(a,i4,a,i4,a)')                  &
            '?line_msg: index ', i,                    &
            ' is out of the Fortran range [1:', n, ']'

       line = errmsg

    else

       line = trim(m%msg(i))

    end if

  end subroutine line_msg


  !
  ! ------------------------------------------------------------------
  ! write_msg
  !
  ! Write the message to a Fortran logical unit.
  ! ------------------------------------------------------------------
  !

  subroutine write_msg(m, n)

    type(submsg), intent(in) :: m
    integer,      intent(in) :: n

    integer :: i

    if (m%num <= 0) then
       return
    end if

    do i = 1, m%num
       write(n, '(a)') trim(m%msg(i))
    end do

  end subroutine write_msg


  !
  ! ------------------------------------------------------------------
  ! print_msg
  !
  ! Print the message to standard output.
  ! ------------------------------------------------------------------
  !

  subroutine print_msg(m)

    type(submsg), intent(in) :: m

    integer :: i

    if (m%num <= 0) then
       return
    end if

    do i = 1, m%num
       print '(a)', trim(m%msg(i))
    end do

  end subroutine print_msg


  !
  ! ------------------------------------------------------------------
  ! get_defmsg
  !
  ! Return pointer to global/default submsg.
  ! ------------------------------------------------------------------
  !

  function get_defmsg() result(m)

    type(submsg), pointer :: m

    if (.not. igot) then

       allocate(default_msg)

       call new_msg(default_msg)

       igot = .true.

    end if

    m => default_msg

  end function get_defmsg


  !
  ! ------------------------------------------------------------------
  ! replace_defmsg
  !
  ! Replace global/default message.
  ! ------------------------------------------------------------------
  !

  subroutine replace_defmsg(m)

    type(submsg), intent(in) :: m

    if (.not. igot) then

       allocate(default_msg)

       call new_msg(default_msg)

       igot = .true.

    end if

    call clear_msg(default_msg)

    call copy_msg(default_msg, m)

  end subroutine replace_defmsg


end module submsg_mod



!
! ======================================================================
!
! C-callable access routines for the default message.
!
! All routines below use bind(C) so ctypes sees the exact symbol names:
!
!     c_clear_defmsg
!     c_print_defmsg
!     c_size_defmsg
!     c_size_line
!     c_resize_defmsg
!     c_trim_defmsg
!     c_line_defmsg
!     c_add_defmsg
!
! ======================================================================
!


!
! ----------------------------------------------------------------------
! c_clear_defmsg
!
! Clear the default error-message buffer.
! ----------------------------------------------------------------------
!

subroutine c_clear_defmsg() bind(C, name="c_clear_defmsg")

  use submsg_mod

  implicit none

  type(submsg), pointer :: d

  d => get_defmsg()

  call clear_msg(d)

end subroutine c_clear_defmsg



!
! ----------------------------------------------------------------------
! c_print_defmsg
!
! Print the default error-message buffer.
! ----------------------------------------------------------------------
!

subroutine c_print_defmsg() bind(C, name="c_print_defmsg")

  use submsg_mod

  implicit none

  type(submsg), pointer :: d

  d => get_defmsg()

  call print_msg(d)

end subroutine c_print_defmsg



!
! ----------------------------------------------------------------------
! c_size_defmsg
!
! Return number of logical lines in the default message.
! ----------------------------------------------------------------------
!

function c_size_defmsg() result(n) bind(C, name="c_size_defmsg")

  use iso_c_binding, only: c_int
  use submsg_mod

  implicit none

  integer(c_int) :: n

  type(submsg), pointer :: d

  d => get_defmsg()

  n = int(size_msg(d), kind=c_int)

end function c_size_defmsg



!
! ----------------------------------------------------------------------
! c_size_line
!
! Return maximum Fortran error-message line length.
! ----------------------------------------------------------------------
!

function c_size_line() result(n) bind(C, name="c_size_line")

  use iso_c_binding, only: c_int
  use submsg_mod

  implicit none

  integer(c_int) :: n

  n = int(nsubmsg, kind=c_int)

end function c_size_line



!
! ----------------------------------------------------------------------
! c_resize_defmsg
!
! ns >= 0:
!     Set storage size.
!
! ns < 0:
!     Ensure abs(ns) free lines.
!
! The integer is intentionally NOT VALUE because existing ctypes callers
! pass a pointer/reference.
! ----------------------------------------------------------------------
!

subroutine c_resize_defmsg(ns) bind(C, name="c_resize_defmsg")

  use iso_c_binding, only: c_int
  use submsg_mod

  implicit none

  integer(c_int), intent(in) :: ns

  type(submsg), pointer :: d

  d => get_defmsg()

  call resize_msg(d, int(ns))

end subroutine c_resize_defmsg



!
! ----------------------------------------------------------------------
! c_trim_defmsg
!
! Remove ns lines from the end of the default message.
!
! The integer is intentionally NOT VALUE because ctypes passes byref().
! ----------------------------------------------------------------------
!

subroutine c_trim_defmsg(ns) bind(C, name="c_trim_defmsg")

  use iso_c_binding, only: c_int
  use submsg_mod

  implicit none

  integer(c_int), intent(in) :: ns

  type(submsg), pointer :: d

  d => get_defmsg()

  call trim_msg(d, int(ns))

end subroutine c_trim_defmsg



!
! ----------------------------------------------------------------------
! c_line_defmsg
!
! Fetch one line from the default message.
!
! i uses C indexing:
!
!       i = 0  -> Fortran line 1
!       i = 1  -> Fortran line 2
!       ...
!
! n is the available C buffer size.
!
! Returned c() is always null terminated when n > 0.
!
! Integers are intentionally NOT VALUE because ctypes passes byref().
! ----------------------------------------------------------------------
!

subroutine c_line_defmsg(i, n, c) bind(C, name="c_line_defmsg")

  use iso_c_binding, only: c_char, c_int, c_null_char
  use submsg_mod

  implicit none

  integer(c_int), intent(in) :: i
  integer(c_int), intent(in) :: n

  character(kind=c_char), intent(out) :: c(*)

  character(len=nsubmsg) :: p

  type(submsg), pointer :: d

  integer :: j
  integer :: nc
  integer :: iline


  !
  ! Nothing can be written to a zero-size buffer.
  !

  if (n <= 0_c_int) then
     return
  end if


  !
  ! Convert C index to Fortran index.
  !

  iline = int(i) + 1


  d => get_defmsg()

  p = ' '

  call line_msg(d, iline, p)


  !
  ! Copy as many characters as fit while leaving room
  ! for the terminating null.
  !

  nc = min(len_trim(p), int(n) - 1)

  if (nc < 0) then
     nc = 0
  end if


  do j = 1, nc
     c(j) = p(j:j)
  end do


  !
  ! Null terminator.
  !

  c(nc + 1) = c_null_char

end subroutine c_line_defmsg



!
! ----------------------------------------------------------------------
! c_add_defmsg
!
! Add a null-terminated C string to the default message.
!
! This performs the C -> Fortran conversion explicitly instead of
! relying on cstring(), avoiding the termination problem encountered
! in the c_lx_open interface.
! ----------------------------------------------------------------------
!

subroutine c_add_defmsg(c) bind(C, name="c_add_defmsg")

  use iso_c_binding, only: c_char, c_null_char
  use submsg_mod

  implicit none

  character(kind=c_char), intent(in) :: c(*)

  character(len=nsubmsg) :: p

  type(submsg), pointer :: d

  integer :: i


  p = ' '


  do i = 1, len(p)

     if (c(i) == c_null_char) then
        exit
     end if

     p(i:i) = c(i)

  end do


  d => get_defmsg()

  call add_msg(d, trim(p))

end subroutine c_add_defmsg
