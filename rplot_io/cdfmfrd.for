      subroutine cdfmfrd(file_id, ifcn, ind, retn)
C
C Read "ifcn" Profile Function into COMMON buffer
C
C 09/24/97  CAL
 
      use datmgr_mod
      use cplotr_mod
      include "netcdf.inc"
 
C Input
      integer file_id
      integer ifcn
      integer ind
C Return
      integer retn
 
C Processing
      integer v_id
      integer ipt, isize
 
      character*21 name
C----------------------------------------------------------------------
C Variable Ids of profile functions are function index +1
C          1st variable = Time
 
      if(lrun_x.eq.0) then
         v_id = ifcn+ncdft
      else
         v_id = ifcn+ncdft_x(lrun_x)
      endif
 
      ipt = locd(ind)
C
      if(lrun_x.eq.0) then
         isize=nzonex(itypr(ifcn)) * ntr
         name=abr(ifcn)
      else
         itype=itypr_x(ifcn,lrun_x)
         izonex=nzonex_x(itype,lrun_x)
         isize=izonex*ntr_x(lrun_x)
         name=abr_x(ifcn,lrun_x)
      endif
C
      nwds(ind) = isize
 
      retn = nf_get_var_real(file_id, v_id, datbuf(ipt))
      if (retn .ne. NF_NOERR) then
         luntrm=lunzer(0)
         write(luntrm,9001) name,NF_STRERROR(retn)
 9001    format(' % CDFMFRD: nf_get_vara_real - ',a/,
     >        '     ',a)
         return
      end if
C
      return
      end
