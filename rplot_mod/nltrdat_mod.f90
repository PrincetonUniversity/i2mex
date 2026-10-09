module nltrdat_mod
  implicit none

  integer, parameter :: Nummax = 3000  ! maximum records in a ...TR.dat file.
  integer, parameter :: Nreclen = 140  ! Length of input records for ...TR.dat
  integer, parameter :: Nreclen2 = Nreclen+2  ! Length of input records for ...TR.dat

  integer :: NNames
  character(len=32), dimension(NumMax) :: Names
  character(len=120), dimension(NumMax) :: Zinput

end module nltrdat_mod
