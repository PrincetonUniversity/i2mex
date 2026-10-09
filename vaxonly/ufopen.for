C******************** START FILE UFOPEN.FOR ; GROUP UFOPEN *************
c
C--------------------------------------------------------------
C  UFOPEN
C  OPEN SEQUENTIAL FILE ON UNIT ILUN
C
      SUBROUTINE UFOPEN(ILUN,ZFILE,ICOD,IER,IBIN)
c
      implicit NONE
c
c
C  ILUN-- UNIT NUMBER
C  ZFILE- CHARACTER FILENAME
C  ICOD-- =1 FOR READ ACCESS
C         =2 FOR WRITE ACCESS
C  IER--- ERROR RETURN (IOSTAT ON OPEN) CODE 0= NORMAL
C
C-------------------------------------------------------------
C  ** modified DMC ** using new portlib routine, is_ascii, to test
C       whether a file is ASCII or not.  I find the "file" command
C       is not to be trusted for portability...
C
C  ** Modified TBT Dec 2000 ** added check for "commands text" as well
C                                          as    " ascii text".
C                              Failing on uscws8.gat.com where file command
C                              returns "commands text".
C                              Added check for "text"
C                              If file protection has "x" turned on,
C                              file command returns "commands text" rather
C                              than "ascii text".
C  ** MODIFIED CAL Aug 1998 ** use "jsystem" don't use csh/isbinary
C
C  ** MODIFIED CAL Jul 1997 ** Check for HDF, netCDF
C
C  ** MODIFIED DMC SEP 1993 ** FOR PORTABILITY, THE OPEN STATEMENTS
C    HAVE BEEN MOVED TO SUBROUTINES GENOPEN AND IOPEN, WHICH ARE
C    CALLED FROM HERE.
C
C  ** MODIFIED DMC APR 1984 **
C    SUPPORT BINARY SEQUENTIAL ACCESS
C
C  argument IBIN -- input if file is OPEN for WRITE (ICOD=2)
C    IBIN-- 0 FOR ASCII
C           1 FOR BINARY ACCESS (could be netCDF or old Ufiles binary or...)
C           2 For a different binary format (Ufiles loss free compressed)
C           3 For HDF
C           4 For NetCDF
C
C  on OPEN for READ (ICOD=1) set IBIN=-1 to have UFOPEN test
C
C              or binary (return IBIN=1)
C        Ufiles new compressed binary (return IBIN=2)
C                  HDF   (return IBIN=3)
C
C
C  on OPEN for READ if IBIN is NOT -1, assume the caller already
C  knows the file is ascii (IBIN=0) or binary (IBIN=1).
C
C  if open for READ and IBIN = -1, the input file type has to
C  be tested.
C
C
      integer ilun,icod,ier,ibin
C
      CHARACTER*(*) ZFILE
C
      integer str_length
C
      integer*1 btest(512), bref(8), hdf_id(4), cdf_id(4)
C
      integer nf_noerr,nf_nowrite
C
      parameter (nf_noerr = 0)
      parameter (nf_nowrite = 0)
C
      integer ilf,i,ios,iascii,idcdf,nf_open,nf_close,irecl
C
      data bref/0,127,0,127,0,127,0,127/
      data hdf_id/   14,     3,   19,    1/
C
C-----------------------------------------------------------------------
C
      ilf=str_length(zfile)
      if(ilf.le.0) then
         write(6,'('' ?UFOPEN:  null filename argument.'')')
         ier=99
         return
      endif
C
C
C  CHECK ARGUMENTS
      IF(ICOD.EQ.1) THEN
C...OPEN for READ
        IF((IBIN.ge.0) .and. (IBIN.le.4)) GO TO 60
        IF(IBIN.EQ.-1) GO TO 10  ! TEST FILE TYPE ASCII OR BINARY
      ELSE IF(ICOD.EQ.2) THEN
C...OPEN for WRITE
         IF((IBIN.ge.0).and.(IBIN.le.4)) GO TO 100
      ENDIF
C
      IER=-1
      WRITE(6,9901) ICOD,IBIN,ZFILE(1:ilf)
 9901   FORMAT(
     >' ?UFOPEN -- INVALID ARGUMENTS, ICOD=',I8,' IBIN=',I8/
     >'  FILE "',A,'" NOT OPENED.')
      RETURN
c-----------------------------------------------------------------------------
C
C-------------------------------
C  TEST FILE TYPE
C
 10   CONTINUE
C
      IBIN=0
C
C ** read access ** determine if file is binary or ascii
C (support binary Ufiles)
C
      call is_ascii(zfile,iascii,ier)
      if(ier.eq.0) then
         if(iascii.eq.1) then
            IBIN=0
         else
            IBIN=1
         endif
      endif
C
C  check for new type of binary Ufile --
C
      if(ibin.eq.1) then
         call cftopen(ilun,zfile(1:ilf),'rb',1,ier)
         if(ier.ne.0) go to 57
         call cftread(ilun,btest,ier)
         if(ier.ne.0) go to 57
         do i=1,8
            if(btest(i).ne.bref(i)) go to 57
         enddo
C
C  OK
         IBIN=2
         go to 58
C
C Check for HDF
 57      do i=1,4
            if(btest(i) .ne. hdf_id(i)) go to 59
         enddo
         IBIN = 3
cxx         type *,'HDF File'
         go to 58
C
C
C     Check for netCDF
 59      continue
         idcdf=-99
         ios = nf_open(zfile(1:ilf),NF_NOWRITE,idcdf)
         if(ios.eq.NF_NOERR) then
            IBIN = 4
            ios=nf_close(idcdf)
cxx            type *,'netCDF File'
         endif
C
 58      continue
         call cftclose(ilun)
      endif
C
C---
C  READ OPEN
C
 60   IOS = 0
C
      IF(IBIN.EQ.0) THEN
        CALL GENOPEN(ILUN,ZFILE,'OLD','ASCII',0,IOS)
      ELSE IF(IBIN.EQ.1) THEN
        CALL GENOPEN(ILUN,ZFILE,'OLD','BINARY',0,IOS)
      ELSE IF(IBIN.EQ.2) THEN
        CALL GENOPEN(ILUN,ZFILE,'OLD','BLOCKD',128,IOS)
      ENDIF
C
      GO TO 500
C
C---
C  WRITE OPEN
C
 100  IOS = 0
C
CCC      write(*,'(a,a)')'In UFOPEN: ZFILE= ',ZFILE
      IF(IBIN.EQ.0) THEN
        CALL GENOPEN(ILUN,ZFILE,'NEW','ASCII',0,IOS)
      ELSE IF(IBIN.EQ.1) THEN
        CALL GENOPEN(ILUN,ZFILE,'NEW','BINARY',0,IOS)
      ELSE IF(IBIN.EQ.2) THEN
        CALL GENOPEN(ILUN,ZFILE,'NEW','BLOCKD',128,IOS)
      ENDIF
C
      GO TO 500
C-------------------------------
 500  CONTINUE
      IER=IOS
      RETURN
      	END
C******************** END FILE UFOPEN.FOR ; GROUP UFOPEN ***************
