C-----------------------------------------------------------------------
C  GENOPEN -- generic OPEN routine
C
C  mod dmc -- support new type, BLOCKD, for sequential access with
C  fixed length records.
C
      subroutine GENOPEN(lun,name,fstat,ftype,irecsz,ios)
C
      integer lun          ! logical unit number
      character*(*) name   ! file name
      character*(*) fstat  ! file status 'OLD' or 'NEW' or 'UNKNOWN'
      character*(*) ftype  ! file type 'ASCII' 'BINARY' or 'DIRECT'
      integer irecsz       ! record size ('DIRECT' files only)
C
C       ** units of irecsz is (4 byte) WORDS on all systems **
C       ** integer function nblkfac converts this to bytes on
C       ** systems for which this is necessary...
C
C output:
      integer ios          ! status code, 0 is normal
C                            ! <0 for argument error (message generated)
C                            ! >0 for file open error
C
C--------------------------------------
C  local declarations
C
      character*10 zstat,ztype
C
C--------------------------------------
C
C  check arguments
C
      IOS=0
C
      zstat=' '
      ztype=' '
      ilstat=len(fstat)
      iltype=len(ftype)
      zstat=fstat(1:ilstat)
      ztype=ftype(1:iltype)
      call uupper(zstat)
      call uupper(ztype)
C
      if( (zstat.ne.'OLD') .and. (zstat.ne.'UNKNOWN') .and.
     >      (zstat.ne.'NEW') ) then
        write(6,'('' ?GENOPEN:  invalid fstat argument:  '',A)') fstat
        IOS=-1
      endif
C
      if( (ztype.ne.'ASCII') .and.
     >	    (ztype.ne.'BINARY') .and.
     >      (ztype.ne.'DIRECT') .and.
     >      (ztype.ne.'BLOCKD') ) then
        write(6,'('' ?GENOPEN:  invalid ftype argument:  '',A)') ftype
        IOS=-2
      endif
C
      if(ztype.eq.'DIRECT') then
        if(irecsz.le.0) then
          IOS=-3
          write(6,9901) irecsz
 9901     format(' ?GENOPEN:  invalid record size for DIRECT file:  ',
     >        i10)
        endif
        iblkfac=nblkfac(ier)
        if(ier.ne.0) then
           write(6,'('' %genopen -- NBLKFAC error on DIRECT access.'')')
           iblkfac=1
        endif
        irecl=iblkfac*irecsz
      endif
C
      if(ios.lt.0) return
C
C  ok, do the file open
C
      if(ztype.eq.'BLOCKD') then
         iln=len_trim(name)
         if(zstat.eq.'OLD') then
            call cftopen(lun,name(1:iln),'rb',1,ios)
         else if(zstat.eq.'NEW') then
            call cftopen(lun,name(1:iln),'wb',1,ios)
         else
            write(6,9902) fstat
 9902       format(
     >         ' ?GENOPEN:  fstat="',a,'" not valid for BLOCKD files.')
         endif
         go to 1000  ! C code open, no fortran error handling.
      endif
C
      ios=iopen(lun,name,zstat,ztype,irecl)
C
 1000 continue
      return
      end
