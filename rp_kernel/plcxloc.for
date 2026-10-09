      subroutine plcxloc(istat,iwrk1,iwrk2,ipx,ztestval,z2,z0)

      use datmgr_mod

C  for the profile currently in the buffer DATBUF(iwrk2...),
C  at each time, find the X array value corresponding to ztestval.
C  if there is more than one matching X value, return z2
C  if there is no matching X value, return z0
C
C---------------------------------
      use cplotr_mod
C---------------------------------
      real zxarr(NR0)
C---------------------------------
C
      inx=nzonex(istat)
C
      IF(.not.NLXFOT(ISTAT)) THEN
         do ix=1,inx
            zxarr(ix)=xarry(ix,istat)
         enddo
      ENDIF
C
      do it=1,ntr
C  X axis at this time
         if(nlxfot(istat)) then
            ipxl=ipx+(it-1)*inx
            do ix=1,inx
               zxarr(ix)=datbuf(ipxl+ix-1)
            enddo
         endif
C  find occurrence(s) of ztestval in the profile
         ipf=iwrk2+(it-1)*inx
         inum=0
         ztes2p=1.0
         do ix=1,inx-1
C  X desired at this time
            zf1=datbuf(ipf+ix-1)
            zf2=datbuf(ipf+ix)
            ztes1=ztestval-zf1
            ztes2=ztestval-zf2
            if((sign(1.0,ztes1)*ztes2).le.0.0) then
               if(zf1.eq.zf2) then
                  inum=inum+2           ! degeneracy
               else
                  if((ztes1.eq.0.0).and.(ztes2p.eq.0.0)) then
                     continue           ! at bdy & already caught
                  else
                     inum=inum+1
                     ixsave=ix
                     zxfac=(ztestval-zf1)/(zf2-zf1)
                     zx1=zxarr(ix)
                     zx2=zxarr(ix+1)
                  endif
               endif
            endif
            ztes2p=ztes2
         enddo
C
C  do interpolation, if there is one and only one location at the test value
C
         if(inum.eq.0) then
            zans=z0
         else if(inum.gt.1) then
            zans=z2
         else
            zans=(1.0-zxfac)*zx1+zxfac*zx2
         endif
C
         datbuf(iwrk1+it-1)=zans
      enddo
C
      return
      end
 
