      subroutine TRDATE(RUNDAT,ORDATE,FILNAM,MODE)
use iso_c_binding, only: fp => c_double
!
!  mod dmc 3 feb 1992 -- standardize the fortran i/o
!
!  --------------------------------------------------------------------------
!
!  MODE=1 :
! ---  Translate date for ordering
!     dd-mmm-yyyy  -->  yyyy-mm-dd   where mm is numerical order of month
!
!  MODE=2 :
! ---  Reverse
!     yyyy-mm-dd  -->  dd-mmm-yyyy
!
!
!  --------------------------------------------------------------------------
!
      character  MONUC(12)*3, MONLC(12)*3, NMON*2, RUNDAT*11, &
                 ORDATE*12, FILNAM*(*)
      DATA  MONUC /'JAN','FEB','MAR','APR','MAY','JUN', &
      		     'JUL','AUG','SEP','OCT','NOV','DEC'/
      DATA  MONLC /'jan','feb','mar','apr','may','jun', &
      		     'jul','aug','sep','oct','nov','dec'/
!
      if(MODE.EQ.1)   THEN
        do I=1,12
          if(   RUNDAT(4:6) .EQ. MONUC(I) .OR. &
                  RUNDAT(4:6) .EQ. MONLC(I))  THEN
            WRITE(NMON,'(I2)',ERR=700)  I
            if(I.LE.9)  NMON='0'//NMON(2:2)
            goto 100
          end if
        end do
        WRITE(6,9001) FILNAM,RUNDAT
 9001   FORMAT(' TRDATE:  File ',A,' date translation error ',A)
        call BAD_EXIT
!
  100   if(RUNDAT(1:1).EQ.' ')  RUNDAT(1:2)='0'//RUNDAT(2:2)
        ORDATE=RUNDAT(8:11)//'-'//NMON//'-'//RUNDAT(1:2)
        return
      else if(MODE.EQ.2) THEN
!***	  DECODE(2,'(I2)',ORDATE(6:7))  I
        Read(ORDATE(6:7),'(I2)')  I     !*** jpj@JET Apr-1995
        RUNDAT=ORDATE(9:10)//'-'//MONUC(I)//'-'//ORDATE(1:4)
        return
      end if
!
      WRITE(6,9002) MODE
 9002 FORMAT(' TRDATE:  bad MODE value',I5)
      call BAD_EXIT
!
!
  700 WRITE(6,'('' TRDATE:  internal write error'')')
      call BAD_EXIT
      end
