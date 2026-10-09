      subroutine cdfnfrd(file_id, ipt, isize, retn)
C
C Read Scalar Functions into COMMON buffer
C
C 09/23/97  CAL
 
      use datmgr_mod
      use cplotr_mod
      include "netcdf.inc"
 
C Input
      integer file_id
      integer ipt           ! pointer into datbuf
C Return
      integer isize         ! size of all Scalar Functions
      integer retn
 
C Processing
      integer v_id
      integer i, ipt1, isize_min
 
C----------------------------------------------------------------------
 
 
C Read Time Dimension  -- this was moved to cdfhrd (dmc 17 Nov 1997)
C====================
C
C Read all Scalar Functions into datbuf
C--------------------------------------
C id of 1st Scalar = NFXT + 2 (includes Time axes # = ncdft)
 
cdbg      type *,'Read Scalar Functions',nft,ntt
 
      luntrm=lunzer(0)
 
      if(lrun_x.eq.0) then
         v_id = nfxt + ncdft
         isize=ntt*nft
         inft=nft
         intt=ntt
      else
         v_id = nfxt_x(lrun_x) + ncdft_x(lrun_x)
         isize=ntt_x(lrun_x)*nft_x(lrun_x)
         inft=nft_x(lrun_x)
         intt=ntt_x(lrun_x)
      endif

      isize_min=20*isize
      if(isize_min.gt.ndbsiz) then
         write(6,*) ' ntt,nft,isize = ',ntt,nft,isize
         write(6,*) ' isize_min,ndbsiz ',isize_min,ndbsiz
         call dmg_datbuf_expand(isize_min)
      endif

      ipt1 = ipt
      do i = 1, inft
         v_id = v_id + 1
         retn = nf_get_var_real(file_id, v_id, datbuf(ipt1))
         if (retn .ne. NF_NOERR) then
            if(lrun_x.eq.0) then
               write(luntrm,9002) abt(i), NF_STRERROR(retn)
            else
               write(luntrm,9002) abt_x(i,lrun_x), NF_STRERROR(retn)
            endif
 9002       format(' % CDFNFRD: nf_get_vara_real - ',a/,
     >           '     ',a)
            return
         end if
         ipt1 = ipt1 + intt
      end do
 
      return
      end
