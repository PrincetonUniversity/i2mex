      subroutine tmklbl(disk,ildisk,dir,ildir,runid,ilrunid,
     >                  zrlbl,ilrlbl)
C
C  dmc 18 Nov 1997 -- trprofil/trscalar support routine
C  get string lengths & make composite label for run
C
C  input:
      character*(*) disk                ! disk
      character*(*) dir                 ! directory (path)
      character*(*) runid               ! TRANSP runid
C
C  output:  the non-blank string lengths (ildisk, etc...), and
C
      character*(*) zrlbl               ! composite label
C
C  local:
C
      character*1 char
C
C-------------------------------
C
C  non-blank string lengths:
C
      ildisk=ilnurd(disk)
      ildir=ilnurd(dir)
      ilrunid=ilnurd(runid)
      if(disk(1:4).eq.'MDS+') then
!
!  MDS+ run access
!
         if(dir.ne.' ') then
            zrlbl='MDS+:'//dir(1:ildir)//'!'//runid(1:ilrunid)
         else
            ildm1=index(disk,'.')-1
            ilamp1=index(disk,'@')+1
            if(min(ildm1,ilamp1).lt.5) then
               zrlbl='MDS+?:'//runid(1:ilrunid)
            else
               zrlbl='MDS+:'//disk(6:ildm1)//'!'//disk(ilamp1:ildisk)//
     >            '!'//runid(1:ilrunid)
            endif
         endif
      else if(dir.eq.' ') then
!
!  just the runid (current directory)
!
         zrlbl=runid
      else
!
!  use tail of dir as part of composite label...
!
         idot=0
         ick=ildir+1
         char=dir(ildir:ildir)
         if(char.eq.'/') ick=ick-1
         if(char.eq.'>') ick=ick-1
         if(char.eq.']') ick=ick-1
         ilast=ick-1
         do ic=1,8
            ick=ick-1
            char=dir(ick:ick)
            if(char.eq.'.') idot=idot+1
            if(idot.eq.2) go to 10
            if(char.eq.' ') go to 10
            if(char.eq.'/') go to 10
            if(char.eq.':') go to 10
         enddo
         ick=ick-1
C
C  fall thru loop is OK
C
 10      continue
         ick=ick+1
         zrlbl=dir(ick:ilast)//'!'//runid(1:ilrunid)
      endif
C
      ilrlbl=ilnurd(zrlbl)
C
      return
      end
