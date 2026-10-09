      subroutine tgethost(zmds,zhost)
c
      use cplotr_mod
c
c  this returns zmds = "MDS+", zhost = <mds+ server machine>
c    --or--
c  zmds = " ", zhost = <current machine> if file based access is
c  being used.
c
      character*(*) zmds               ! path (directory) argument
      character*(*) zhost              ! runid argument
c
c--------------------------------------
      integer indx
c
      if((lfdisk.gt.0).and.(fdisk(1:4).eq.'MDS+')) then
         zmds='MDS+'
         indx=index(fdisk,'@')-1
         if(indx.le.5) then
            zhost=' ?unknown '
         else
            zhost=fdisk(6:indx)
         endif
      else
c
c  file access
c
         zmds=' '
         zhost=' '
         call hostnm(zhost)          ! gets current host
c
      endif
c
      return
      end
