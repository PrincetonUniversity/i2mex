C******************** START FILE DMIDEL.FOR ; GROUP DATMGR ******************
C--------------------------------------------------------
C  DMIDEL
C
C  DELETE A DATA ENTRY J AND MAKE A FREE SLOT
C
      SUBROUTINE DMIDEL(J)
C
C	Last updated:
C	     6/25/91  TBT  Added check of J>0.
C 
      use datmgr_mod
C
      IF (J .LE. 0)  RETURN    ! Nothing to be done - Undefined.
C
      if(no_delete) then
         call errmsg_exit(
     >      ' ?? DMIDEL: DATBUF entry deletions prohibited!')
      endif
C
      LAVAIL=MIN0(LAVAIL,J)
      NWDS(J)=0
      LOCD(J)=0
      MPRIO(J)=0
C  PATCH CHAIN, BOTH WAYS
      LNEXT(LPREV(J))=LNEXT(J)
      LPREV(LNEXT(J))=LPREV(J)

#ifdef __DEBUG
      call dmprin('dmidel',1000)
#endif

      RETURN
      END
C******************** END FILE DMIDEL.FOR ; GROUP DATMGR ******************
