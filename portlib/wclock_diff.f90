subroutine wclock_diff(ivals_start, ivals_end, t_seconds_r4)

  ! single precision version -- compute elapsed wall clock time, from
  ! two time stamps, in seconds

  ! for full comments see the R8 precision version, below...

  implicit NONE

  !--------------------
  ! INPUT:

  integer, dimension(8), intent(in) :: ivals_start  ! start time
  !  (returned from f90 date_and_time call)

  integer, dimension(8), intent(in) :: ivals_end    ! stop time
  !  (returned from f90 date_and_time call)

  real, intent(out) :: t_seconds_r4
  !  elapsed wall clock time (seconds)

  !------------------------
  INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
  real(kind=R8) :: t_seconds_r8
  !------------------------

  call wclock_diff_r8(ivals_start,ivals_end,t_seconds_r8)
  t_seconds_r4 = t_seconds_r8

end subroutine wclock_diff

subroutine wclock_diff_r8(ivals_start, ivals_end, t_seconds)
  
  ! subroutine to difference integer data from two calls to f90 intrinsic 
  ! routine DATE_AND_TIME:
  !   integer, dimension(8) :: ivals
  !   call DATE_AND_TIME(values=ivals)  On Sept. 25, 2007 at 9:26:15.794 am
  !      ivals(1) = year, e.g.   2007
  !      ivals(2) = month, e.g.  9
  !      ivals(3) = day, e.g.    25
  !      ivals(4) = (time zone adjustment, ignored)
  !      ivals(5) = hour, e.g.   9
  !      ivals(6) = minute, e.g. 26
  !      ivals(7) = second, e.g. 15
  !      ivals(8) = msecs, e.g.  794

  ! contained routine wclock_diff returns the elapsed time between two
  ! calls to DATE_AND_TIME, the data captured in separate 8 word integer
  ! arrays.  The elapsed time can be returned in units specified by the
  ! user:  seconds, minutes, hours, days, years (assuming 365.25 days/year)
  !        or months (years*12.0).

  ! ACCURACY: the module assumes a leap year once every 4 years; millenial
  ! corrections are ignored (a later version may take this into account).

  ! NOTE also: for now, ivals(4) (time zone adjustment) is IGNORED; it is
  ! assumed without checking that the time zone adjustment is unchanged 
  ! in the two time stamps being compared.

  ! DMC 25 Sep 2007
  !-------------------------------

  implicit NONE

  INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)

  real*8, parameter :: dpyav = 365.25_R8 ! days/year (approx.)

  real*8, parameter :: x1000 = 1000.0_R8 ! msecs/sec
  real*8, parameter :: x1000th = 0.001_R8 ! secs/msec
  real*8, parameter :: x60   = 60.0_R8   ! minutes/hour, seconds/minute
  real*8, parameter :: x24   = 24.0_R8   ! hours/day
  real*8, parameter :: x12   = 12.0_R8   ! months/year

  real*8, parameter :: x100  = 100.0_R8
  real*8, parameter :: ONE   = 1.0_R8
  real*8, parameter :: ZERO  = 0.0_R8

  ! normal year days/month

  integer, dimension(12), parameter :: dpmon_norm = &
       (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

  ! leap year days/month

  integer, dimension(12), parameter :: dpmon_leap = &
       (/ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

  !--------------------
  ! INPUT:

  integer, dimension(8), intent(in) :: ivals_start  ! start time
  !  (returned from f90 date_and_time call)

  integer, dimension(8), intent(in) :: ivals_end    ! stop time
  !  (returned from f90 date_and_time call)

  !--------------------
  ! OUTPUT:  

  real(kind=R8), intent(out) :: t_seconds  ! elapsed time in seconds

  !--------------------
  ! LOCAL:
  integer ia1(8),ia2(8),isign,iyear
  real(kind=R8) :: zh1,zh2
  real(kind=r8) :: zs1,zs2
  real(kind=r8) :: tdiff,daysum,tot_day1,day1,day2
  real(kind=r8) :: tdiff_hours

  !--------------------
  !  first see if ivals_start is BEFORE or AFTER 

  isign = 0

  do
     zh1 = x1000*x1000*ivals_start(1) + x100*x100*ivals_start(2) + &
          x100*ivals_start(3) + ONE*ivals_start(5)

     zh2 = x1000*x1000*ivals_end(1) + x100*x100*ivals_end(2) + &
          x100*ivals_end(3) + ONE*ivals_end(5)

     if(zh1.gt.zh2) then
        isign = -1
        ia1 = ivals_end
        ia2 = ivals_start

     else if(zh1.lt.zh2) then
        isign = 1
        ia1 = ivals_start
        ia2 = ivals_end

     else
        zs1= x60*ivals_start(6) + ONE*ivals_start(7) + ivals_start(8)/x1000
        zs2= x60*ivals_end(6) + ONE*ivals_end(7) + ivals_end(8)/x1000

        if(zs1.gt.zs2) then
           isign = -1
           ia1 = ivals_end
           ia2 = ivals_start

        else if(zs1.lt.zs2) then
           isign = 1
           ia1 = ivals_start
           ia2 = ivals_end

        else
           isign = 1
           tdiff = ZERO
           exit ! the two times are identical...
        endif
     endif

     !  OK the two times are in ia1 & ia2
     !  ia1 guaranteed earlier than ia2.

     day1 = days_in(ia1) ! days into the year ia1(1)
     day2 = days_in(ia2) ! days into the year ia2(2)

     if(ia1(1).eq.ia2(1)) then
        !  two time stamps in SAME year
        tdiff = day2-day1
        exit
     endif

     daysum = ZERO

     !  first count days in complete years between the two times...

     do iyear = ia1(1)+1,ia2(1)-1
        if(leap(iyear)) then
           daysum = daysum + ONE*366
        else
           daysum = daysum + ONE*365
        endif
     enddo
     
     !  number of days in year containing 1st time stamp
     if(leap(ia1(1))) then
        tot_day1 = ONE*366
     else
        tot_day1 = ONE*365
     endif

     !  two time stamps in DIFFERENT years
     tdiff = daysum + (tot_day1 - day1) + day2
     exit
  enddo

  !-------------------------
 
  t_seconds = isign*tdiff*x24*x60*x60

  contains
    ! PRIVATE routines
    
    !--------------------------------
    logical function leap(iyear)

      !  return T if leap year (every 4th year)

      !----------------
      integer, intent(in) :: iyear
      !----------------

      !----------------
      integer :: iy4
      !----------------

      iy4 = iyear/4
      iy4 = iy4*4

      leap = (iy4.eq.iyear)

    end function leap

    !--------------------------------
    real(kind=R8) function days_in(ia8)
      
      !  return the time into the year ia8(1) in time stamp ia8(1:8), in days

      !----------------
      integer, intent(in) :: ia8(8)
      !----------------

      !----------------
      real(kind=R8) :: dsum
      integer, dimension(12) :: montab
      integer :: imon
      !----------------

      if(leap(ia8(1))) then
         montab = dpmon_leap
      else
         montab = dpmon_norm
      endif

      ! count full months prior to month of time stamp

      dsum = ZERO
      do imon = 1,ia8(2)-1
         dsum = dsum + ONE*montab(imon)
      enddo

      ! add in the number of days of month of time stamp

      dsum = dsum + ONE*ia8(3)

      ! +hours

      dsum = dsum + ia8(5)/x24

      ! +minutes

      dsum = dsum + ia8(6)/(x24*x60)

      ! +seconds + msecs

      dsum = dsum + (ONE*ia8(7) + x1000th*ia8(8))/(x24*x60*x60)

      !-------------

      days_in = dsum

    end function days_in

end subroutine wclock_diff_r8
