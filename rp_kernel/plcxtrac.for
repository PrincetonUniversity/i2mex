      subroutine plcxtrac(istat,iwrk1,iwrk2,ipx,imin,ixof)
C
C  extract one of:  max(f),min(f),x(max(f)),x(min(f)) from profile function
C  f at datbuf(iwrk2...); store result at datbuf(iwrk1...)
C
C  **note** this produces a function of time on the **profile** timebase;
C  use "ttintrp" to map back to the **scalar** functions timebase.
C
      use datmgr_mod
      use cplotr_mod
C
      integer istat                     ! profile function type
      integer iwrk1                     ! DATBUF address:  where to write
      integer iwrk2                     ! DATBUF address:  where f is
      integer ipx                       ! DATBUF address:  where x is
C
      integer imin                      ! =1:  want min, not max
      integer ixof                      ! =1:  want x, not f
C
      logical ilxfot
C
C----------------------------------------------
C
      if(istat.le.0) then
        call zermsg(' ?plcxtrac:  no profile function, istat.le.0')
      else
        inx=nzonex(istat)
        ilxfot=nlxfot(istat)
C
        do it=1,ntr
           ia=iwrk1+it-1
           iaf=iwrk2+(it-1)*inx
           zflim=datbuf(iaf)
           if(ilxfot) then
              iax=ipx+(it-1)*inx
              zxlim=datbuf(iax)
           else
              zxlim=xarry(1,istat)
           endif
           do ix=2,inx
              zf=datbuf(iaf+ix-1)
              if(ilxfot) then
                 zx=datbuf(iax+ix-1)
              else
                 zx=xarry(ix,istat)
              endif
              if(imin.eq.1) then
                 if(zf.lt.zflim) then
                    zflim=zf
                    zxlim=zx
                 endif
              else
                 if(zf.gt.zflim) then
                    zflim=zf
                    zxlim=zx
                 endif
              endif
           enddo
C
           if(ixof.eq.1) then
              datbuf(ia)=zxlim
           else
              datbuf(ia)=zflim
           endif
C
        enddo
      endif
C
      return
      end
