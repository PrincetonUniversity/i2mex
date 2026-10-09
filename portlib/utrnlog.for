      subroutine UTRNLOG(table,loglnam,physnam,lphys,istat)
c
c  translate logical name to physical name
c
c  arguments:
      character*(*) table	!table to search
      character*(*) loglnam	!logical name to be translated
      character*(*) physnam	!(return) physical name
      integer       lphys	!(return) length of physical name
      integer       istat	!(return) status code
C
C  ISTAT=1 DENOTES SUCCESS
C
C  on UNIX system:
C   table is ignored
C   loglnam is the environment variable name
C   physnam is the value returned
C   lphys is the length w/o trailing blanks of the value returned
C   istat is 1 unless lphys=0, in which case istat is 2 on exit.
C
      ist=1
      if(loglnam(1:1).eq.'$') ist=2
      call mpi_sget_env(loglnam(ist:),physnam,iertmp)
      lphys=len_trim(physnam)
      if(lphys.eq.0) then
        istat=2
      else
        istat=1
      endif
C
      return
C
      end
