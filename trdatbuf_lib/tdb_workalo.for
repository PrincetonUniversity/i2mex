C---------------------------------------------------------
C  TDB_WORKALO
C
C  ALLOCATE SPACE IN WORK DATA BUFFER
C
C  BOMB OUT IF NONE IS AVAILABLE
C
C  INPUT ISIZE = AMOUNT OF SPACE (NO. WORDS) NEEDED
C  OUTPUT ILOC = LOCATION OF FIRST WORD OF ALLOCATED SPACE
C
C  COMMON QUANTITY LFREE_W IS INCREMENTED TO NEXT FREE SEGMENT
C
      SUBROUTINE TDB_WORKALO(d,ISIZE,ILOC)
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
      IEND=D%LFREE_W+ISIZE-1
C
C  check if need to expand buffer
C
      call workbuf_expand(d,iend,isizb)
C
C  ALLOCATE SPACE in buffer
C
      ILOC=D%LFREE_W
      D%LFREE_W=D%LFREE_W+ISIZE
C
      RETURN
      END
