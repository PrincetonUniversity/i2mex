SUBROUTINE find_io_unit(ilun)

  ! Return next available Fortran I/O Unit Specifier

  ! unit number returned in range MIN_ILUN (120) to MAX_ILUN (999)
  ! return ilun = -1 if none are available but it is probably not
  ! possible for ~1000 files to be open simultaneously for other system 
  ! reasons... (note: MIN_ILUN increased to 120 as pathscale compiler
  ! reserves units 100:102 and perhaps others will as well).

  IMPLICIT none
  INTEGER, INTENT(OUT) :: ilun
  LOGICAL :: open
  INTEGER, PARAMETER :: MIN_ILUN=120, MAX_ILUN=999

  integer :: iu,ifound

  ilun = -1
  ifound = 0

  DO iu=MIN_ILUN, MAX_ILUN
     INQUIRE(UNIT=iu, OPENED=open)
     IF (.NOT. open) THEN
        ifound=1
        EXIT
     ENDIF
  END DO

  if(ifound.eq.1) ilun = iu

END SUBROUTINE find_io_unit

SUBROUTINE find_io_unit_list(n,ilun_list)

  ! Return next N available Fortran I/O Unit Specifiers

  ! unit number returned in range MIN_ILUN (120) to MAX_ILUN (999)
  ! return ilun_list(j)=-1 if none are available for #j and up.

  IMPLICIT none
  INTEGER, INTENT(IN)  :: n  ! number of available LUNs desired.
  INTEGER, INTENT(OUT) :: ilun_list(n)
  LOGICAL :: open
  INTEGER, PARAMETER :: MIN_ILUN=120, MAX_ILUN=999

  integer :: iu,ifound

  ilun_list(1:n) = -1
  ifound = 0

  DO iu=MIN_ILUN, MAX_ILUN
     INQUIRE(UNIT=iu, OPENED=open)
     IF (.NOT. open) THEN
        ifound=ifound+1
        ilun_list(ifound) = iu
        if(ifound.eq.n) exit
     ENDIF
  END DO

END SUBROUTINE find_io_unit_list
