      subroutine plcdfdt(jtyp,iwrk1,iwrk2,ipt)

      use datmgr_mod

C
C  dmc Aug 1999 -- calculate time derivative
C
C  datbuf(iwrk1...) = workspace -- overwritten, copy of orig. data
C  datbuf(iwrk2...) = the data to be differentiated
C                     data is replaced with time derivative of data
C
C  jtyp -- the data x axis type (-1 if a scalar fcn of time)
C
C  ipt  -- ptr to scalar workspace -- needed for eldot correction of
C          time derivative of flux surface functions
C
C----------------------------
C
      use cplotr_mod
C
C--------------------------------------
C
      if(jtyp.eq.0) then
         call zermsg(' ??plcdfdt -- jtyp=0 on call.')
         call abortt
      endif
C
      if(jtyp.eq.-1) then
         inx=1
         int=ntt
      else
         inx=nzonex(jtyp)
         int=ntr
      endif
C
      isiz=inx*int
      call copyr4(datbuf(iwrk2),datbuf(iwrk1),isiz)
C
      do it=1,int
C
         ia1=iwrk1+(it-1)*inx
         ia2=iwrk2+(it-1)*inx
C
         itm=max(1,it-1)
         itp=min(int,it+1)
         iam=iwrk1+(itm-1)*inx
         iap=iwrk1+(itp-1)*inx
C
         call idtcalc(datbuf(ia1),datbuf(iap),datbuf(iam),
     >        datbuf(ia2),inx,1.0,it,itp,itm,IPT,jtyp)
      enddo
C
      return
      end
 
