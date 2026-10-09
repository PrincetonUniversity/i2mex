C---------------------------------------------------------
C  TDB_DATALO
C
C  ALLOCATE SPACE IN PHYSICS (MEASURED) DATA BUFFER
C  (TRANSP COMMON BLOCK 5.2)
C
C  BOMB OUT IF NONE IS AVAILABLE
C
C  INPUT ISIZE = AMOUNT OF SPACE (NO. WORDS) NEEDED
C  OUTPUT ILOC = LOCATION OF FIRST WORD OF ALLOCATED SPACE
C
C  COMMON QUANTITY LFREE IS INCREMENTED TO NEXT FREE SEGMENT
C
      SUBROUTINE TDB_DATALO(d,ISIZE,ILOC)
C
      use trdatbuf_obj
      IMPLICIT NONE
      type (trdatbuf) :: d          ! data buffer
      integer, intent(in) :: isize  ! size of chunk to allocate
      integer, intent(out) :: iloc  ! address of chunk allocated
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER iend,isizb
!============

! 12/04/08 CLF: don't allocate size=0 
      if (isize .le. 0) return
!
      IEND=d%LFREE+ISIZE
C
C  check if need to expand buffer
C
      call datbuf_expand(d,iend,isizb)
C
C  ALLOCATE SPACE in buffer
C
      ILOC=d%LFREE
      d%LFREE=d%LFREE+ISIZE
C
      RETURN
      END
