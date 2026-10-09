      subroutine rp_dmgalo(isize,ind,iprio)
c
c  RPLOT datmgr space allocation interface routine.
c    * always allocate with priority =5
c    * assigned passed priority if successful
c
      use datmgr_mod
c
      call dmgalo(isize,ind,5)
      mprio(ind)=iprio
c
      return
      end
