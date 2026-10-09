!-----------------------------------------------------------------------
!  PARSID -- PARSE OUT TRANSP RUN ID AND TOKAMAK FROM COMBINED INPUT
!
! History:
! 10/16/96 CAL : Support RunId with 4 - 6 digit shot number
!                                   2 - 3 digit run  number
!
      subroutine PARSID(STRING,NRUN,TOK)
use iso_c_binding, only: fp => c_double
!
      character(len=*) STRING	! INPUT STRING E.G. "12345A01PBXM"
!				! COULD BE OLD STYLE "1234TFTR"
!                               ! e.g. 123456A123TFTR
      character(len=*) NRUN	! OUTPUT RUN ID, "12345A01" OR "1234"
      character(len=*)  TOK	! OUTPUT TOKAMAKD ID, "PBXM" OR "TFTR"
!
      character(len=10) :: ZDIGIT
!
!  TWO RUN ID "STYLES" ARE SUPPORTED, THE OLD 4 DIGIT RUN NUMBER STYLE
!  AND THE NEW 8 character nnnnnAmm SHOT-TRY NOMENCLATURE
!
!------------------------------------
!
      TOK='????'
!
      if(len(nrun).lt.10) then
         ilen=len(nrun)
         write(6,9901) ilen
 9901    format( &
             ' ?parsid:  len(rundid)=',i2, &
             ' incompatible with 6 digit shot number.')
         nrun='len_error'
         return
      end if
!
      ZDIGIT='0123456789'
!
      NRUN='ERROR!!!'
!
      ILEN=INDEX(STRING,' ')-1
      if(ILEN.LE.0) ILEN=LEN(STRING)
!
!  TEST FOR OLD VS. NEW STYLE IS WHETHER THE 5TH character OF THE INPUT
!  STRING IS A DIGIT:  DIGIT INDICATES NEW STYLE
!
      if(ILEN.LT.5) goto 1000
!
      IDIG=INDEX(ZDIGIT,STRING(5:5))
!
      if(IDIG.EQ.0) THEN
!
!  OLD STYLE
!
!  CHECK LENGTH
!
        if((ILEN.LT.7).OR.(ILEN.GT.8)) goto 1000
!
!  OK, CHECK THAT 1ST FOUR CHARS ARE DIGITS
!
        do 10 IC=1,4
          if(INDEX(ZDIGIT,STRING(IC:IC)).EQ.0) goto 1000
 10     continue
!
!  GOOD
!
        NRUN=STRING(1:4)
        TOK=STRING(5:ILEN)
!
      else
!
!  NEW STYLE
!
!  CHECK valid Shot Number
!
        ic = 0
        do 20 while (idig .ne. 0)
          ic = ic+1
          IDIG = INDEX(ZDIGIT,STRING(IC:IC))
 20     continue
        i1 = ic - 1
        if (i1 .lt. 5 .or. i1 .gt. 6) goto 1000
          if (i1 .eq. 6 .and. string(1:1) .eq. '0') goto 1000

!
! Check valid Run Number
          idig = 1
        do 21 while (idig .ne. 0)
          ic = ic+1
          IDIG = INDEX(ZDIGIT,STRING(IC:IC))
 21     continue
        i2 = ic - i1 - 2
        if (i2 .lt. 2 .or. i2 .gt. 3)  goto 1000

        NRUN=STRING(1:ic-1)
        TOK=STRING(ic:ILEN)
!
      end if
!
 1000 continue
      return
      end
