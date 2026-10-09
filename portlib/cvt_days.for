      subroutine cvt_days(iseconds,iarray)
      integer iseconds                  ! number of seconds, input (.gt.0)
      integer iarray(4)                 ! equiv. no. of days, hours, mins, secs
C
C  on output,  iarray(1) = no. of complete 24 hour days covered by iseconds
C              iarray(2) = no. of hours (past a complete day) covered
C                          btw 0 and 23
C              iarray(3) = no. of minutes (past a complete hour) covered
C                          btw 0 and 59
C              iarray(4) = no. of seconds (past a cmoplete minute) covered
C                          btw 0 and 59
C
      data isecmin/60/                  ! seconds / minute
      data isechour/3600/               ! seconds / hour
      data isecday/86400/               ! seconds / day
C
      iwork=iseconds
C
      iarray(1)=iwork/isecday
      iwork=iwork-isecday*iarray(1)
C
      iarray(2)=iwork/isechour
      iwork=iwork-isechour*iarray(2)
C
      iarray(3)=iwork/isecmin
      iarray(4)=iwork-isecmin*iarray(3)
C
      return
      end
