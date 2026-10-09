      subroutine mfinq(MFILN,NZONES,IBLKSZ,IRECSZ,MFBLKI)
C
C input:
      character*(*) MFILN               ! filename
      integer NZONES                    ! old format recordsize
      integer IBLKSZ                    ! new format recordsize
C
C output:
      integer IRECSZ                    ! actual recordsize
      logical MFBLKI                    ! .TRUE. for new format
C
C  broken out of rplot_sub/inirun.for -- dmc -- 17 Nov 1997
C  inquire as to format of MF.PLN file
C
C============
      IRECSZ=0
 
C============
C CAL 07/30/98
C Note:
C On UNIX systems, INQUIRE does not return RECL
C
      ibytfac=nblkfac(ier)
      IBLKL=ibytfac*IBLKSZ
      IBLKOLD=ibytfac*NZONES
C
      IF((IRECSZ.EQ.IBLKL).OR.(IRECSZ.NE.IBLKOLD)) THEN
C
C  ASSUME 512 BYTE NEW FORMAT MF FILE
         MFBLKI=.TRUE.
         IRECSZ=IBLKL
      ELSE
C
C  OLD FORMAT MF FILE
         IRECSZ=IBLKOLD
         MFBLKI=.FALSE.
      ENDIF
C
C  dmc -- divide ibytfac back out
C
      irecsz=irecsz/ibytfac
C
      return
      end
