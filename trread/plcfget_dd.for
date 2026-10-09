      subroutine plcfget_dd(zpath,zdisk,zdir,ier)
c
c  make disk/dir args for trprofil/trscalar
c
      character*(*) zpath
      character*(*) zdisk
      character*(*) zdir
      integer ier
C
      character*5 ztest5
C
      ztest5=zpath(1:min(len(zpath),5))
      call uupper(ztest5)
      if(ztest5.eq.'MDS+:') then
         call plcfget_mds(zpath,zdisk,zdir,ier)
         if(ier.ne.0) return
      else
         zdisk=' '
         zdir=zpath
      endif
      return
      end
