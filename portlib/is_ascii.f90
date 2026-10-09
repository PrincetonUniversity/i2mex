subroutine is_ascii(filename,iascii,ier)
  use iso_c_binding, only: c_char
 
  !  determine if a file is ascii
 
  character*(*), intent(in) :: filename
  integer, intent(out) :: iascii          ! =1: ascii, =0: no, =-1: error
  integer, intent(out) :: ier             ! =0: OK, otherwise "C" errno code.
 
  ! --------------------------
 
 
  character(kind=c_char) cfname(1+len(filename))
 
  ! --------------------------
 
  call cstring(filename(1:len_trim(filename)),cfname,'2C')
 
  call portlib_isascii(cfname,iascii,ier)   ! call the "C" routine to do it...
 
  return
 
end subroutine is_ascii
