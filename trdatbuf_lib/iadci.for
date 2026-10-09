C  IADCI
C  Devel. from IADC3 -- replace ntmmx with NTIME2
C   Calc. the pointer to an element in the time-interpolated
C                                      FRBC moments C3 array
C       JFS  25 May 95
C
      Function iadci(d, ibase, ita, ixa, ima, iia )
 
      use trdatbuf_obj
      IMPLICIT NONE
      type (trdatbuf) :: d
 
C---------------------------------------------------------------------
 
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ibase,ita,ixa,ima,iia,iadci,in3,isizb
!============
      in3   = (ima + 1) + (iia-1) * (d%mmax+1)
 
      iadci = ibase + (ita-1)
     +              + (ixa-1) * D%NTIME2
     +              + (in3-1) * D%NTIME2 * d%nxmmx
 
      call datbuf_expand(d,iadci,isizb)
 
      Return
 
      End
! 19jan2003 fgtok -s r8_precision.sub all.sub "r8con.csh conversion"
