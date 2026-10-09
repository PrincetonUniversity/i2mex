C******************** START FILE PLMCALC.FOR ; GROUP PLOTMM ******************
C-----------------------------------------------------------------------
C  PLMCALC
C
C  CALCULATE sin and cos tables for use in PLMGEO, PLMAGEO, PLMMCV
C                                                & PLAMMCV
 
 
      SUBROUTINE PLMCALC(Jmom, Itheta)
 
 
C  INPUT:
C        Jmom   = # of moments (not counting 0)
C        Itheta = # of theta values to calculate.
 
      use cplotr_mod
      use plfmpa_mod
 
      real ztheta(NaxMMP)
 
c---------------------------
 
      if(Itheta.gt.NaxMMP) then
         write(lunzer(0),9901) Itheta
 9901    format(' ?PLMCALC:  Itheta = ',i10/
     1	  '  exceeds CPLOTR COMMON parameter NaxMMP')
         ncostabl=0
         return
      else
         ncostabl=Itheta
      endif
 
 
      ZTwoPi = 6.2831853072
 
      DO 20 ITH=1,Itheta
 
        ZTHETA(ith)=ZtwoPi*(ITH-1)/FLOAT(Itheta-1)
 
 20   CONTINUE
C
      call plmcalc_tbl(NaxMom,Jmom,itheta,Ztheta,ZCosTabl,ZSinTabl)
C
      RETURN
      END
c-----------------------
      subroutine plmcalc_lim(imaxth)
c
c  return max no. of theta points
c
      use cplotr_mod
      use plfmpa_mod
c
      imaxth=NaxMMP
c
      return
      end
c-----------------------
      subroutine plmcalc_tbl(mxmom,jmom,itheta,ztheta,zcostabl,zsintabl)
c
c  compute sin-cos table for equilibrium lookup
c  mxmom - max no. of moments in table;
c  jmom  - actual no. of momehts
c  itheta - no. of theta points
c
      real ztheta(itheta)               ! theta grid
      real zcostabl(0:mxmom,itheta)     ! cos(im*theta) table
      real zsintabl(0:mxmom,itheta)     ! sin(im*theta) tableb
c
      do ith = 1,itheta
         zth=ztheta(ith)
         DO IM=0,Jmom
 
            ZCosTabl(Im,Ith) = COS(IM*ZTH)
            ZSinTabl(Im,Ith) = SIN(IM*ZTH)
 
         enddo
      enddo
c
      return
      end
C******************** END FILE PLMCALC.FOR ; GROUP PLOTMM ******************
