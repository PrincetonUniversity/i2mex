c--------------------------------------------------------------
c  get the equilibrium all the data points at one time
c
      subroutine pltammt(jtyp,ztime,zdelta,zbufr,zbufz,id1,id2)

      use datmgr_mod
      use cplotr_mod
      use plfmpa_mod

c  input:
      integer jtyp                      ! x-axis id code
      real ztime                        ! time at which equilibrium needed
c
      integer id1,id2                   ! array dimensions
c
c  output:
      real zbufr(id1,id2),zbufz(id1,id2) ! the equilibrium
c
c   zbufr(naxmmp,inx),zbufz(naxmmp,inx)
c   R(theta,x),Z(theta,x)
c-------------------
C
C  LOCAL--
C
      REAL ZRSMOM(0:NAXMOM),ZYSMOM(0:NAXMOM)      ! Sin    moments
      REAL ZRCMOM(0:NAXMOM),ZYCMOM(0:NAXMOM)      ! Cosine moments
C
      real zztime                       ! local time
      integer iexp                      ! range flag
      integer ix                        ! surface index, no. of surfaces
      integer im                        ! moment index
      integer i1,i2,ipx1,ipx2           ! datbuf addresses
c
c-------------------
C  NON-TEMPORAL AXIS
      INX=NZONEX(JTYP)
C
      if((id1.gt.naxmmp).or.(id2.ne.inx)) then
         call zermsg(' %pltammt:  unexpected array dimensions.')
         do i2=1,id2
            do i1=1,id1
               zbufr(i1,i2)=0.0
               zbufz(i1,i2)=0.0
            enddo
         enddo
         return
      endif
c
      call pltmmfx(ztime,zdelta,jtyp)
c
      ztime1=ztime-zdelta
      ztime2=ztime+zdelta
c
      if((zdelta.le.0.0).or.(ztime+zdelta.le.time3(1)).or.
     >     (ztime-zdelta.ge.time3(ntr)) .or. ntr<2 .or. 
     >     ztime1>=ztime2) then
c
c  interpolate at one time
c
         zztime=ztime
         CALL PLTIMI(ZZTIME,IT1,IT2,Z1,Z2,IEXP)
C
C  GET THE MOMENTS FOR EACH X VALUE
C
         DO IX=1,INX
            DO IM=0,IMOM
C  R MOMENTS
               I1=NDPTR(ILCAMM(1,IM),INX,IT1)+IX-1
               I2=NDPTR(ILCAMM(1,IM),INX,IT2)+IX-1
               ZRCMOM(IM)=DATBUF(I1)*Z1+DATBUF(I2)*Z2 ! Cosine
 
               I1=NDPTR(ILCAMM(2,IM),INX,IT1)+IX-1
               I2=NDPTR(ILCAMM(2,IM),INX,IT2)+IX-1
               ZRSMOM(IM)=DATBUF(I1)*Z1+DATBUF(I2)*Z2 ! Sin
C  Y MOMENTS
               I1=NDPTR(ILCAMM(3,IM),INX,IT1)+IX-1
               I2=NDPTR(ILCAMM(3,IM),INX,IT2)+IX-1
               ZYCMOM(IM)=DATBUF(I1)*Z1+DATBUF(I2)*Z2 ! Cosine
 
               I1=NDPTR(ILCAMM(4,IM),INX,IT1)+IX-1
               I2=NDPTR(ILCAMM(4,IM),INX,IT2)+IX-1
               ZYSMOM(IM)=DATBUF(I1)*Z1+DATBUF(I2)*Z2 ! Sin
            enddo
C
C  DEFINE THE CONTOUR
            CALL PLAMMCV(ZRCMOM,ZRSMOM,ZYCMOM,ZYSMOM,IMOM,IX,
     >         zbufr,zbufz,id1,id2)
C
         enddo
      else
c
c  time average
c
         do ix=1,inx
            do im=0,imom
               zrcmom(im)=0.0
               zrsmom(im)=0.0
               zycmom(im)=0.0
               zysmom(im)=0.0
            enddo
c
            zwsum=0.0
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
                  do im=0,imom
C  R MOMENTS
                     I1=NDPTR(ILCAMM(1,IM),INX,IT1)+IX-1
                     I2=NDPTR(ILCAMM(1,IM),INX,IT2)+IX-1
                     ZRCMOM(IM)=zrcmom(im)+zdtw*
     >                  (DATBUF(I1)*Z1+DATBUF(I2)*Z2) ! Cosine
 
                     I1=NDPTR(ILCAMM(2,IM),INX,IT1)+IX-1
                     I2=NDPTR(ILCAMM(2,IM),INX,IT2)+IX-1
                     ZRSMOM(IM)=zrsmom(im)+zdtw*
     >                  (DATBUF(I1)*Z1+DATBUF(I2)*Z2) ! Sin
C  Y MOMENTS
                     I1=NDPTR(ILCAMM(3,IM),INX,IT1)+IX-1
                     I2=NDPTR(ILCAMM(3,IM),INX,IT2)+IX-1
                     ZYCMOM(IM)=zycmom(im)+zdtw*
     >                  (DATBUF(I1)*Z1+DATBUF(I2)*Z2) ! Cosine
 
                     I1=NDPTR(ILCAMM(4,IM),INX,IT1)+IX-1
                     I2=NDPTR(ILCAMM(4,IM),INX,IT2)+IX-1
                     ZYSMOM(IM)=zysmom(im)+zdtw*
     >                  (DATBUF(I1)*Z1+DATBUF(I2)*Z2) ! Sin
                  enddo
               endif
            enddo
c
            do im=0,imom
               zrsmom(im)=zrsmom(im)/zwsum
               zrcmom(im)=zrcmom(im)/zwsum
               zysmom(im)=zysmom(im)/zwsum
               zycmom(im)=zycmom(im)/zwsum
            enddo
C
C  DEFINE THE CONTOUR
            CALL PLAMMCV(ZRCMOM,ZRSMOM,ZYCMOM,ZYSMOM,IMOM,IX,
     >         zbufr,zbufz,id1,id2)
C
         enddo
      endif
C
      return
      end
