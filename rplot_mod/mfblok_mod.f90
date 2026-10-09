!---------------------------------------------------------------------
!  MFBLOK - RPLOT MF FILE READ COMMON BLOCK
!
!  dmc -- multiple run support added 14 Nov 1997
!
!  ** CPLOTR ** must also be declared, before MFBLOK
!
!      use cplotr_mod
!      include 'MFBLOK'
!
module mfblok_mod
  use cplotr_mod, only: naxxtra
  implicit none

  integer, parameter :: IBLKSZ = 128
  !
  LOGICAL MFBLKI			! block format flag
  INTEGER MFLUNI			! LUN
  INTEGER,allocatable,dimension(:) :: MFHDR ! header data
  REAL MFDATA(IBLKSZ),MFDUM(IBLKSZ)	! data buffers
  !
  logical MFBLKI_X(naxxtra)	! block format flags
  integer MFLUNI_X(naxxtra)	! LUNs
  integer ,allocatable,dimension(:,:) :: MFHDR_X ! headers

end module mfblok_mod
