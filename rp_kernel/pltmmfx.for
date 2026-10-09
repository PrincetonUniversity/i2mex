      subroutine pltmmfx(ztime,zdelta,jtyp)

      use datmgr_mod
      use cplotr_mod
      use plfmpa_mod

c  get the (time averaged or time interpolated) x axis, store in
c  array XF in common
c
c  input:
c
      real ztime                        ! time of interest, seconds
      real zdelta                       ! +/- zdelta (seconds)
      integer jtyp                      ! x axis id code
c
c  if zdelta.le.0 -- interpolate to ztime
c  if zdelta.gt.0 -- average over ztime-zdelta to ztime+zdelta
c
C
C  LOCAL--
C
      real zztime                       ! local copy of time var.
      integer iexp                      ! interpolation info
      integer it,ix                     ! loop indices
      integer ipx1,ipx2                 ! datbuf addreses
C
c-------------------------------------------------------------
      IF(.NOT.NLXFOT(JTYP)) THEN
c
c  fixed x axis
c
         do ix=1,nzonex(jtyp)
            XF(IX)=XARRY(IX,JTYP)
         enddo
c
         return
c
      endif
c
c  time varying x axis
c
      inxi=nzonex(jtyp)
c
      if((zdelta.le.0.0).or.(ztime+zdelta.le.time3(1)).or.
     >   (ztime-zdelta.ge.time3(ntr))) then
c
c  interpolation (or @ endpt)
c
         zztime=ztime
         CALL PLTIMI(ZZTIME,IT1,IT2,Z1,Z2,IEXP)
C  NON-TEMPORAL AXIS
         DO IX=1,INXI
            IPX1=NXPTR(JTYP,IT1)+IX-1
            IPX2=NXPTR(JTYP,IT2)+IX-1
            XF(IX)=DATBUF(IPX1)*Z1+DATBUF(IPX2)*Z2
         ENDDO
c
      else
c
c  time average
c
         ztime1=ztime-zdelta
         ztime2=ztime+zdelta
         if(.not.nlxfot(jtyp)) then
            do ix=1,inxi
               xf(ix)=xarry(ix,jtyp)
            enddo
         else
            zwsum=0.0
            do ix=1,inxi
               xf(ix)=0.0
            enddo
            do it=1,ntr-1
               zta=time3(it)
               ztb=time3(it+1)
               ztest1=max(ztime1,zta)
               ztest2=min(ztime2,ztb)
               if(ztest2.gt.ztest1) then
                  zztime=0.5*(ztest1+ztest2)
                  zdtw=ztest2-ztest1
                  zwsum=zwsum+zdtw
                  CALL PLTIMI(ZZTIME,IT1,IT2,Z1,Z2,IEXP)
                  do ix=1,inx
                     IPX1=NXPTR(JTYP,IT1)+IX-1
                     IPX2=NXPTR(JTYP,IT2)+IX-1
                     XF(IX)=xf(ix)+zdtw*
     >                  (DATBUF(IPX1)*Z1+DATBUF(IPX2)*Z2)
                  enddo
               endif
            enddo
            do ix=1,inxi
               xf(ix)=xf(ix)/zwsum
            enddo
         endif
      endif
C
      return
      end
