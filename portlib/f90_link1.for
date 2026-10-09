	subroutine f90_link1
C
C  this is a token f90 routine.  In mixed f77/f90 environments, sometimes
C  an actual f90 routine must be linked, in order for the load process to
C  go correctly.
C
C  in portlib this is called out of "err_end".
C
      integer ivals(8)
C
      character*8 zzdate
      character*10 ztime
      character*5 zone
C
      call date_and_time(zzdate,ztime,zone,ivals)
      if(ivals(2).gt.12) write(6,1001) ivals(2)
 1001	format(' f90_link1:  ?? month .gt. 12:  ',i10)
C
      return
      end
