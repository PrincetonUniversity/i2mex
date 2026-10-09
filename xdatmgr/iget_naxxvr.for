      subroutine iget_naxxvr(ival,inr0)

      ! internal routine for xdatmgr: return NAXXVR value from 'CPLOTR'
      ! INCLUDE file (COMMONs); also return NR0

      use cplotr_mod

      ival = NAXXVR
      inr0 = NR0

      return
      end
