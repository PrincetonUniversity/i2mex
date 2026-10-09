      subroutine mgdelfun(ipkg,ifun,ipos)
C
      use cplotr_mod
C
      integer ipkg                      ! package #
      integer ifun                      ! function #
      integer ipos                      ! position w/in mg
C
C  (previously part of muldef) .. delete function #ifun from package #ipkg
C  error checking is responsibility of caller!
C
C  if the last function is deleted the whole package is deleted!
C
      NFUNS=INFB(IPKG)
      INFB(IPKG)=INFB(IPKG)-1
      NFUNS=NFUNS-1
      itypb=iintb(ipkg)
C
      lunt=lunzer(0)
      IF((ITYPB.EQ.0))
     >     WRITE(lunt,2042) IFUN,ABR(IFUN),abb(ipkg)
      IF((ITYPB.EQ.1)) WRITE(lunt,2042) IFUN,ABT(IFUN),abb(ipkg)
 2042 FORMAT(' FUNCTION # ',I3,', "',A,'" DELETED FROM PACKAGE "',a,'"')
C
      IF(IPOS.GT.NFUNS) GO TO 1000
      DO 428 IK=IPOS,NFUNS
        IKP1=IK+1
        IFUNB(IK,IPKG)=IFUNB(IKP1,IPKG)
 428  CONTINUE
C
 1000 continue
      if(nfuns.eq.0) then
         write(lunt,
     >      '('' %mgdelfun:  multigraph package empty, removed:  '',a)')
     >      abb(ipkg)
         call mgdelpkg(ipkg)
      endif
      return
      end
