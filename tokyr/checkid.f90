!-----------------------------------------------------------------------
!  CHECKID -- PARSE OUT TRANSP RUN ID AND TOKAMAK FROM COMBINED INPUT
!
! History:
! 10/16/96 CAL : Support RunId with 4 - 6 digit shot number
!                                   2 - 3 digit run  number
! 01/07/97 CAL : Fix lrid
! 01/06/98 CAL : Fix extracting shot number
!
subroutine checkid(zstring, lrid, lshot)
  implicit none
  ! Argument variables:
  character(len=*), intent(in) :: zstring ! INPUT STRING E.G. "12345A01PBXM"
                                          ! COULD BE OLD STYLE "1234TFTR"
                                          ! e.g. 123456A123TFTR
  integer, intent(out) :: lrid, lshot     ! return values: length of RunId/Shot
  ! Local variables:
  character(len=10) :: zdigit
  integer :: ilen,idig,ic,i1,i2

  ! TWO RUN ID "STYLES" ARE SUPPORTED, THE OLD 4 DIGIT RUN NUMBER STYLE
  ! AND THE NEW 8 CHARACTER nnnnnAmm SHOT-TRY NOMENCLATURE

  lshot = 0
  zdigit = '0123456789'
  !
  ilen = index(zstring,' ')-1
  IF(ilen.LE.0) ilen=len(zstring)
  !	type 1,  string(1:ilen), ilen
  ! 1	format(1x,a20,i3)
  if (ilen .eq. 4 ) go to 1200
  if (ilen .eq. 5 ) go to 1300
  !
  !  TEST FOR OLD VS. NEW STYLE IS WHETHER THE 5TH CHARACTER OF THE INPUT
  !  STRING IS A DIGIT:  DIGIT INDICATES NEW STYLE
  !
  !  CHECK valid shot number
  !
  idig = 1
  ic = 0
  do while (idig .ne. 0 .and. ic .lt. ilen)
    ic = ic+1
    idig = index(zdigit,zstring(IC:IC))
  end do
  i1 = ic - 1
  !if (i1 .eq. 4 .and. i1 .eq. ilen) go to 1100
  if (ic .eq. 4 ) go to 1200
  if (i1 .eq. 4 ) go to 1100
  if (i1 .lt. 5 .or. i1 .gt. 6) go to 1000
  if (i1 .eq. 6 .and. zstring(1:1) .eq. '0') go to 1000
  !
  ! Check valid run number
  !
  idig = 1
  do while (idig .ne. 0 .and. ic .lt. ilen)
    ic = ic+1
    idig = index(zdigit,zstring(IC:IC))
  end do
  if (ic .eq. ilen .and. idig .ne. 0) ic = ic+1
  i2 = ic - i1 - 2
  if (i2 .lt. 2 .or. i2 .gt. 3)  go to 1000
  !
1100 lrid = ic - 1
  lshot = i1
  return

1200 lrid = 4
  lshot = 4
  return

1300 lrid = 5
  lshot = 5
  return

1000 lrid = 0
  return
end subroutine checkid
