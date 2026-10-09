      subroutine plcxntrp(iwrk1,iwrk2,istat,ipf,zconst,zxcept)
C
C  interpolate profile data in datbuf(iwrk2...) to x data in time
C  dependent scalar function at datbuf(ipf...) (if ipf.gt.0) or to
C  x=const (if ipf.eq.0).
C
C  if x is out of range, the resulting value is zxcept.
C  write result to datbuf(iwrk1...)
C
C---------------------------------
      use datmgr_mod
      use cplotr_mod
C---------------------------------
      real zxarr(NR0)
C---------------------------------
C
      inx=nzonex(istat)
C
      IF(NLXFOT(ISTAT)) THEN
         CALL DMGXOT(ISTAT,IND1,IND2)
         IPX=LOCD(IND1)
         IF(ISTAT.EQ.2) IPX=LOCD(IND2)
      else
         do ix=1,inx
            zxarr(ix)=xarry(ix,istat)
         enddo
      ENDIF
C
      itt=1
      do it=1,ntr
         zt=time3(it)
C  X axis at this time
         if(nlxfot(istat)) then
            ipxl=ipx+(it-1)*inx
            do ix=1,inx
               zxarr(ix)=datbuf(ipxl+ix-1)
            enddo
         endif
C  X desired at this time
         if(ipf.eq.0) then
            zx=zconst
         else
            if(zt.le.time(1)) then
               zx=datbuf(ipf)
            else if(zt.ge.time(ntt)) then
               zx=datbuf(ipf+ntt-1)
            else
C  find index time zone
 10            continue
               if(zt.gt.time(itt+1)) then
                  itt=itt+1
                  go to 10
               endif
C  interpolate to profile time from scalar time
               if(time(itt).eq.time(itt+1)) then
                  zfac1=0.5             ! degenerate case
               else
                  zfac1=(zt-time(itt))/(time(itt+1)-time(itt))
               endif
               zfac1=max(0.0,zfac1)
               zfac0=1.0-zfac1
               zx=datbuf(ipf+itt-1)*zfac0+datbuf(ipf+itt)*zfac1
            endif
         endif
C
C  OK -- find interpolation in x space.  CAUTION:  x might not be monotonic.
C
         ixcep=0
         zdelx2p=1.0
C
         zdxmin=abs(zx-zxarr(1))
         idxmin=1
         zfxmin=0.0
C
         zxmin=zxarr(inx)
         zxmax=zxarr(inx)
C
         do ix=1,inx-1
C
            zxmin=min(zxmin,zxarr(ix))
            zxmax=max(zxmax,zxarr(ix))
C
            zdelx1=zx-zxarr(ix)
            zdelx2=zx-zxarr(ix+1)
C
            if(abs(zdelx2).lt.zdxmin) then
               zdxmin=abs(zdelx2)
               idxmin=ix
               zfxmin=1.0
            endif
            if((sign(1.0,zdelx1)*zdelx2).le.0.0) then
C  zx is btw the two zxarr values
               if(zxarr(ix).eq.zxarr(ix+1)) then
                  ixcep=ixcep+2         ! degeneracy (seems unlikely)
               else
                  if((zdelx1.eq.0.0).and.(zdelx2p.eq.0.0)) then
                     continue           ! at bdy & caught already
                  else
                     ixcep=ixcep+1
                     ixsave=ix
                     zxfac=(zx-zxarr(ix))/(zxarr(ix+1)-zxarr(ix))
                  endif
               endif
            endif
            zdelx2p=zdelx2
         enddo
C
C  do interpolation, if there is one and only one location at the test value
C
         if(ixcep.gt.1) then
            zans=zxcept                 ! exception (non-unique)
         else if(ixcep.eq.0) then
            if(zdxmin.lt.1.e-6*max(abs(zxmin),abs(zxmax))) then
               ixcep=1
               ixsave=idxmin
               zxfac=zfxmin
            else
               zans=zxcept              ! exception (test val not found)
            endif
         endif
C
         if(ixcep.eq.1) then
            zf1=datbuf(iwrk2+inx*(it-1)+ixsave-1)
            zf2=datbuf(iwrk2+inx*(it-1)+ixsave)
            zans=(1.0-zxfac)*zf1 + zxfac*zf2
         endif
C
         datbuf(iwrk1+it-1)=zans
      enddo
C
      return
      end
 
