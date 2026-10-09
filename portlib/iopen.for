C-----------------------------------------------------------------------
C  IOPEN -- collection of OPEN statements
C   this routine should not be called directly -- call via subroutine
C   GENOPEN -- the calling arguments are explained there.
C     argument checking is done in GENOPEN.
C
      integer function IOPEN(lun,fname,fstat,ftype,irecsz)
C
      integer lun          ! logical unit number
      character*(*) fname  ! file name
      character*(*) fstat  ! file status 'OLD' or 'NEW'
      character*(*) ftype  ! file type 'ASCII' 'BINARY' or 'DIRECT'
      integer irecsz       ! record size ('DIRECT' files only)
C
      ISTAT=0
C
C           the IOSTAT=ISTAT expression is used in the OPEN
C    statements to capture the error status of the OPEN operation, and
C    this is returned as the function value to the caller.
C
      inb=index(fname,' ')-1
      if(inb.le.0) inb=len(fname)
C
C
      if(ftype.eq.'ASCII') then
C
        if(fstat.eq.'OLD') then
C
          OPEN(unit=lun,file=fname(1:inb),status='OLD',action='READ',
     >            IOSTAT=ISTAT)
C
        else
C
          OPEN(unit=lun,file=fname(1:inb),status=fstat,
     >            access='SEQUENTIAL',IOSTAT=ISTAT)
C
        endif
C
      else if(ftype.eq.'BINARY') then
C
        if(fstat.eq.'OLD') then
C
          OPEN(unit=lun,file=fname(1:inb),status='OLD',action='READ',
     >            access='SEQUENTIAL',form='UNFORMATTED',
     >            IOSTAT=ISTAT)
C
        else
C
          OPEN(unit=lun,file=fname(1:inb),status=fstat,
     >            access='SEQUENTIAL',form='UNFORMATTED',
     >            IOSTAT=ISTAT)
C
        endif
C
      else if(ftype.eq.'DIRECT') then
C
        if(fstat.eq.'OLD') then
C
          OPEN(unit=lun,file=fname(1:inb),status='OLD',action='READ',
     >            access='DIRECT',recl=irecsz,
     >            IOSTAT=ISTAT)
C
        else
C
          OPEN(unit=lun,file=fname(1:inb),status=fstat,
     >            access='DIRECT',recl=irecsz,
     >            IOSTAT=ISTAT)
C
        endif
C
C
      else if(ftype.eq.'BLOCKD') then
C
        if(fstat.eq.'OLD') then
C
          OPEN(unit=lun,file=fname(1:inb),status='OLD',action='READ',
     >            access='SEQUENTIAL',recl=irecsz,
     >            form='UNFORMATTED',IOSTAT=ISTAT)
C
        else
C
          OPEN(unit=lun,file=fname(1:inb),status=fstat,
     >            access='SEQUENTIAL',recl=irecsz,
     >            form='UNFORMATTED',IOSTAT=ISTAT)
C
        endif
C
      endif
C
      IOPEN=ISTAT
C
      return
      end
