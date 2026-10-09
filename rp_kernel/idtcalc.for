C******************** START FILE IDTCALC.FOR ; GROUP IDTCALC ******************
C------------------------------
C  IDTCALC
 
      SUBROUTINE IDTCALC(ZF1,ZFP,ZFM,ZF2,INX,ZSGN,IT,ITP,ITM,
     >    IPT,JTYP)
 
C  created dmc Feb 1996 -- based on IXCALC
C
C  evaluate time derivative (at fixed flux) of flux surface function
C   incorporate eldot*grad(f) correction
C
C  d/dt(f)|flux = d/dt(f)|xi - xi*eldot*df/dxi
C
C  CAPABILITY TO USE TIME VARYING GEOMETRY INFO FROM DATA BUFFER HAS
C  BEEN ADDED.  SEE SUBROUTINE DMGGEO (PLDMGR.FOR)
C
C  INPUT VARIABLES:
C
C  ZF1  LOCAL ARRAY CONTAINING FUNCTION at current time
C  ZFP  LOCAL ARRAY containing function at next time
C  ZFM  LOCAL ARRAY containing function at prev time
C  INX  DIMENSIONALITY OF ZF
C
C  ZSGN - =1.0 OR -1.0:  MULTIPLY PROFILE BY THIS
C
C  IT= INDEX TO CURRENT TIME-- USED FOR TIME DEPENDENT GEOMETRY
C  ITP,ITM -- indices at next and preceding time, for finite diff eval.
C
C  IPT= ptr to scalar data (to find ELDOT data if needed)
C  JTYP=1 x axis type of function; 1=zone ctr'd, 2=bdy ctr'd
C
C  OUTPUT:  ZF2-- time differentiated TRANSFORMED PROFILE
 
      use datmgr_mod

      use cplotr_mod
C
      DIMENSION ZF1(INX),ZF2(INX)
      dimension zfm(inx),zfp(inx)
C
      dimension zdravfac(NR0)  ! gradient correction dmc Feb 96
      real zcorr(NR0)
      real zxi(NR0)
C
      LOGICAL MOMRUN
 
      Real ZLmax
      Data ZLmax /1.e15/
 
      DATA it0/1/
C
C-----------------------------------------------------------------------
C
C  IF GEOMETRY IS TIME DEPENDENT, evaluate d/dt correction (ELDOT factor)
C
      IF(NLTGEO.and.(LELDOT.GT.0).and.(LX.GT.0).and.(LXB.GT.0)
     >     .and.(NTT.gt.0).and.((jtyp.eq.1).or.(jtyp.eq.2))) THEN
C  look up ELDOT at current time
         ztime=time3(it)
         it0p=it0+1
         if(ztime.lt.time(1)) then
            it0=1
            it0p=2
            zfac=0.0
         else if(ztime.gt.time(ntt)) then
            it0=ntt-1
            it0p=ntt
            zfac=1.0
         else
 20         continue
            it0p=it0+1
            if(ztime.lt.time(it0)) then
      	 it0=it0-1
      	 go to 20
            endif
            if(ztime.gt.time(it0p)) then
      	 it0=it0+1
      	 go to 20
            endif
C
            zfac=(ztime-time(it0))/(time(it0p)-time(it0))
C
         endif
         ipeldot=ipt+(LELDOT-1)*NTT
         ipeld1=ipeldot+it0-1
         ipeld2=ipeld1+1
C
         zeldot=datbuf(ipeld1)+zfac*(datbuf(ipeld2)-datbuf(ipeld1))
C
C  look up XI values
C
         if(JTYP.EQ.1) then
            ip=ndptr(lx,inx,it)
         else if(JTYP.eq.2) then
            ip=ndptr(lxb,inx,it)
         endif
         call copyr4(datbuf(ip),zxi,inx)
C
C  compute zcorr = -xi*eldot* df/dXI
C
         inxm1=inx-1
         do ix=2,inxm1
            ixm=ix-1
            ixp=ix+1
            zdfdxi=(zf1(ixp)-zf1(ixm))/(zxi(ixp)-zxi(ixm))
            zcorr(ix)=-zxi(ix)*zeldot*zdfdxi
         enddo
C  near the axis...
         if(jtyp.eq.1) then
            zfac=0.5
         else
            zfac=0.6666667
         endif
         zdfdxi=zfac*(zf1(2)-zf1(1))/(zxi(2)-zxi(1))
         zcorr(1)=-zxi(1)*zeldot*zdfdxi
C  near the edge...
         zdfdxi=(zf1(inx)-zf1(inxm1))/(zxi(inx)-zxi(inxm1))
         zcorr(inx)=-zxi(inx)*zeldot*zdfdxi
C
      else
         do ix=1,inx
            zcorr(ix)=0.0
         enddo
      endif
C
C----------------------------------
C
C  main time derivative
C
      if(itm.eq.itp) then
         do ix=1,inx
            zf2(ix)=0.0
         enddo
      else
         if(jtyp.gt.0) then
            zdti=1.0/(time3(itp)-time3(itm))
         else
            zdti=1.0/(time(itp)-time(itm))
         endif
         do ix=1,inx
            zf2(ix)=zdti*(zfp(ix)-zfm(ix)) + zcorr(ix)
            zf2(ix)=zsgn*zf2(ix)
         enddo
      endif
C
C  THATS ALL FOLKS
C
      RETURN
      END
C******************** END FILE IDTCALC.FOR ; GROUP IDTCALC ******************
