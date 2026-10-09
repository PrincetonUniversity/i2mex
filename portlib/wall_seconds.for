!-----------------------------------------------------------------------
! return the elapsed time, in seconds, since 00:00, Jan. 1, 2000.
! (the algorithm is copied from Charles Karney's "wall_time()"
!
      SUBROUTINE wall_seconds(wall_secs)
      integer wall_secs        ! wall clock seconds since 1/1/2000
!
      integer time(8)
      integer y,m,d
!
!   Charles' description:
! Return wall clock time as the number of seconds since 00:00, January 1,
! 2000.  Use the standard rule for the Gregorian calendar---leap years every
! multiple of 4, unless it's a multiple of 100 and not a muliple of 400.
! Arithmetic is simplified by moving January and February to the previous
! year.  The pattern of month lengths from March thru January is regular with
! a 5-month period---31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31.  The 5-month
! period is 153 days long.  On a computer with 32-bit integers this won't
! overflow until 2068.
!------------------------------------
!
      isecs=0
!
      call clocaltim(time)
!
      y=time(1)+int((time(2)+9)/12)-1 ! Move January and February to prev. year
      m=mod(time(2)+9,12)             ! March is now month 0
!     CAL: .for code is processed as fixed format
      d=365*y+int(y/4)-int(y/100)+int(y/400)+int((m+4)*153/5)+time(3)
     >  -730548
!           January 1, 2000 is day 0
 
      wall_secs=60*(60*(24*d+time(5))+time(6))+time(7)
!
      return
      end
