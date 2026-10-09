module tr_getnl
!
  save
!
!  namelist text
!
  integer :: nltext_status = 99  ! -1 -- not read, 0 -- read, 1 -- error
  integer :: nltext_nlines = 0   ! #of lines in namelist text
  integer :: mltext_nlines = 0   ! #of lines in expanded namelist text
  character(120), dimension(:), allocatable :: nltext
  character(120), dimension(:), allocatable :: mltext
  integer, dimension(:), allocatable :: nltext_lens
  integer, dimension(:), allocatable :: mltext_lens
!
end module
