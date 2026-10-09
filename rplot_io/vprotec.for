C-------------------------------------------------------------------
      subroutine vprotec(tokstr,istat)
C
C  vprotec -- verify access permission
C
C    dmc Mar 1995 -- except at PPPL this routine should not do anything!
C      here at PPPL it is requested to use this to prevent general
C      TRANSP user access to JT60 data.
C
C    rewritten DMC Nov 1996 -- access list now in a file.
C      assume wide open if cannot find file.
C
      implicit NONE
C
C-------------------------------------
C
      character*(*) tokstr  ! tok or tok.yy string, input
      integer istat         ! status code, output. istat=0:  OK
C
C-------------------------------------
C
      integer numacl,icln,il,itok,ilnb
C
      character*80 zbuff,zline
C
      character*4 tokacl(100)
      logical     grantd(100)
C
      integer str_length,ilen
C
      data numacl/-1/
C
      save tokacl
      save grantd
C
C-------------------------------------
C
      istat=0
C
      return
      end
