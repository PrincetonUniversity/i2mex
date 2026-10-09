      subroutine pltmmt(jtyp,ztime,zdelta,zbufr,zbufz,inth,inxi)

      use datmgr_mod
      use cplotr_mod
      use plfmpa_mod

C  return equilibrium in arrays R(theta,x), Z(theta,x)
C
C  input:
      integer jtyp                      ! xb axis code
      real ztime                        ! target time
      real zdelta
      integer inth,inxi                 ! R,Z array dimensions
C
C  output:
      real zbufr(inth,inxi),zbufz(inth,inxi)
C
C  LOCAL--
C
      REAL ZRMOM(NAXMOM),ZYMOM(NAXMOM)  ! moments arrays @ ztime
C
      integer iexp                      ! extrapolation flag
C
      integer ix,i1,i2                  ! loop indices
      integer ipx1,ipx2                 ! datbuf addresses
C
C----------------------------------------------------
C
      if((inth.gt.naxmmp).or.(inxi.ne.nzonex(jtyp))) then
         call zermsg(' %pltmmt:  unexpected array dimensions.')
         do i2=1,inxi
            do i1=1,inth
               zbufr(i1,i2)=0.0
               zbufz(i1,i2)=0.0
            enddo
         enddo
         return
      endif
c
c  get x axis
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
c  interpolation (or @ endpt)
c
         zztime=ztime
         CALL PLTIMI(ZZTIME,IT1,IT2,Z1,Z2,IEXP)
C
C  GET THE MOMENTS FOR EACH X VALUE
C
         DO IX=1,INXI
            I1=NDPTR(ILR0,INXI,IT1)+IX-1
            I2=NDPTR(ILR0,INXI,IT2)+IX-1
            ZRMOM0=DATBUF(I1)*Z1+DATBUF(I2)*Z2
C
            DO IM=1,IMOM
C  R MOMENTS
               I1=NDPTR(ILCMM(1,IM),INXI,IT1)+IX-1
               I2=NDPTR(ILCMM(1,IM),INXI,IT2)+IX-1
               ZRMOM(IM)=DATBUF(I1)*Z1+DATBUF(I2)*Z2
C  Y MOMENTS
               I1=NDPTR(ILCMM(2,IM),INXI,IT1)+IX-1
               I2=NDPTR(ILCMM(2,IM),INXI,IT2)+IX-1
               ZYMOM(IM)=DATBUF(I1)*Z1+DATBUF(I2)*Z2
            enddo
C  DEFINE THE CONTOUR
            CALL PLMMCV(ZRMOM0,ZRMOM,ZYMOM,IMOM,IX,
     >         zbufr,zbufz,inth,inxi)
C
         enddo
      else
c
c  time average
c
         do ix=1,inxi
            zrmom0=0.0
            do im=1,imom
               zrmom(im)=0.0
               zymom(im)=0.0
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
                  I1=NDPTR(ILR0,INXI,IT1)+IX-1
                  I2=NDPTR(ILR0,INXI,IT2)+IX-1
                  ZRMOM0=zrmom0+zdtw*(DATBUF(I1)*Z1+DATBUF(I2)*Z2)
C
                  DO IM=1,IMOM
C     R MOMENTS
                     I1=NDPTR(ILCMM(1,IM),INXI,IT1)+IX-1
                     I2=NDPTR(ILCMM(1,IM),INXI,IT2)+IX-1
                     ZRMOM(IM)=zrmom(im)+zdtw*
     >                  (DATBUF(I1)*Z1+DATBUF(I2)*Z2)
C     Y MOMENTS
                     I1=NDPTR(ILCMM(2,IM),INXI,IT1)+IX-1
                     I2=NDPTR(ILCMM(2,IM),INXI,IT2)+IX-1
                     ZYMOM(IM)=zymom(im)+zdtw*
     >                  (DATBUF(I1)*Z1+DATBUF(I2)*Z2)
                  enddo
 
               endif
            enddo
c
            zrmom0=zrmom0/zwsum
            do im=1,imom
               zrmom(im)=zrmom(im)/zwsum
               zymom(im)=zymom(im)/zwsum
            enddo
c
C  DEFINE THE CONTOUR
            CALL PLMMCV(ZRMOM0,ZRMOM,ZYMOM,IMOM,IX,
     >         zbufr,zbufz,inth,inxi)
c
         enddo
c
      endif
C
      return
      end
