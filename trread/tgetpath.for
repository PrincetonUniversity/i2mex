      subroutine tgetpath(zpath,zrunid)

      use cplotr_mod
c
c  this returns zpath = "MDS+", zrunid = " "
c    --or--
c  zpath=(disk:)directory, zrunid = runid
c
c  if the former, the run tree is open and ancillary information
c  e.g. the TRANSP namelist can be read directly;
c
c  if the latter, zpath & zrunid give the info necessary to form the
c  namelist filename:
c     <runid>TR.DAT if zpath is blank
c     <zpath><runid>TR.DAT if zpath is not blank
c
      character*(*) zpath               ! path (directory) argument
      character*(*) zrunid              ! runid argument
c
c--------------------------------------
c
      if((lfdisk.gt.0).and.(fdisk(1:4).eq.'MDS+')) then
         zpath='MDS+'
         zrunid=' '
      else
c
c  file access
c
         zrunid=runid
         zpath=' '
         if (lfdir .gt. 0) zpath=fdir(1:lfdir)
c
      endif
c
      return
      end
