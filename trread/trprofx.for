      subroutine trprofx(zdisk, zdir, zrunid, zabbrev, zxname,
     >   inxnew,intimes,ztime,zfdat,zxdat,ier)
c
c  read a pair of functions f & x, for which space has been allocated...
c  use trprofil
c
c  arguments as in trprofil...
c
      implicit NONE
c
      character*(*) zdisk
      character*(*) zdir
      character*(*) zrunid
      character*(*) zabbrev
      character*(*) zxname              ! for 2nd fcn
c
      integer inxnew,intimes            ! dimensions...
      real ztime(intimes)
      real zfdat(inxnew,intimes)
      real zxdat(inxnew,intimes)
c
      integer ier
c
c local:
c
      integer isize,ier0,idum
      character*64 zlabel
      character*32 zunits
      integer itype
      integer inxdum,intdum
      integer lunzer
c
c--------------------------------------------
c
      ier0=ier
      isize=inxnew*intimes
C
      ier=-99                           ! keep tree open
      call trprofil(zdisk, zdir, zrunid, zabbrev, INTIMES, ISIZE,
     1             zlabel, zunits, itype, inxdum,
     2             intdum, ztime, zfdat, ier)
C
C  exit now on error
C
      if(ier.ne.0) then
         write(lunzer(0),9902) ier
 9902    format(' %trprofx:  trprofil returned error code:',i5)
         call tconnect_close(idum)
         go to 1000
      endif
C
C  now get x axis data
C
      ier=ier0                          ! close tree when done (or not)
      call trprofil(zdisk, zdir, zrunid, zxname, INTIMES, ISIZE,
     1             zlabel, zunits, itype, inxdum,
     2             intdum, ztime, zxdat, ier)
C
C  exit now on error
C
      if(ier.ne.0) then
         write(lunzer(0),9903) ier
 9903    format(' %trprofx:  trprofil returned error code:',i5)
         call tconnect_close(idum)
         go to 1000
      endif
C
      CALL TIMCK1(ZTIME,INTIMES)
C
 1000 continue
      return
      end
 
