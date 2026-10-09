C******************** START FILE DMGXOT.FOR ; GROUP PLDMGR ******************
C---------------------------------------------------------------
C  DMGXOT
C
C  READ IN TIME-VARYING X AXIS INTO DATA AREA, OR, UPDATE ACCESS
C  CODE IF ALREADY IN DATA AREA
C
      SUBROUTINE DMGXOT(JTYPX,IND1,IND2)
C
      use datmgr_mod
      use cplotr_mod
C
C  JTYPX-- X AXIS ID CODE (INPUT)
C
C  IND1,IND2-- X AXIS DATA LOCATIONS (OUTPUT)
C
C  IF AXIS TYPE IS ZONE BDY/CTR, THE COMPLIMENTARY CTR/BDY DATA IS
C  ALSO SOUGHT
C
C  EXIT IF DATA IS NOT A FCN OF TIME
C
      IND1=0
      IND2=0
      IF(NLTRANSP .AND. .NOT.NLXFOT(JTYPX)) RETURN
C
C  READ AUXILLIARY AXIS (NOT ZONE BDY OR CTR COORDINATE)
C
      IF(NLTRANSP .AND. JTYPX.LT.3) GO TO 10
C
      CALL DMDLOC(ABR(NFX(JTYPX)),IND1,ISIZ,IPT)
      IF(IND1<=0) THEN
         CALL DMGFXT(NFX(JTYPX),IND1)
         CALL DMGXOT_SPREDM(JTYPX,IND1)
      END IF
      IF (.NOT. NLTRANSP) IND2=IND1   ! needed for cdfcon logic
C
      RETURN
C
C  READ ZONE AXES -- this is a TRANSP run with time varying x axes...
C
 10   CONTINUE
C
C  first read in and copy X and XB if necessary...
C
      CALL DMDLOC('%XXXC',IND1,ISIZ,IPT)
      if(ind1.eq.0) then
         ifcn=0
         do i=1,nfxt
            if(abr(i).eq.'X') then
               ifcn=i
               exit
            endif
         enddo
         CALL DMGFXT(ifcn,IND1)       ! "X"
         iloc=locd(ind1)
         inum=nwds(ind1)
         CALL DMDLOC('%X_XC',IND1,ISIZ,IPT)
         ilocx=locd(ind1)
         inumx=nwds(ind1)
         if(inum.ne.inumx) then
            call errmsg_exit(
     >         '? ZC X(%XXXC) axis size inconsistency in dmgxot!')
         endif
         datbuf(ilocx:ilocx+inumx-1)=datbuf(iloc:iloc+inum-1)
         DMGLBL(IND1)='%XXXC'
      endif
      CALL DMDLOC('%XXXB',IND2,ISIZ,IPT)
      if(ind2.eq.0) then
         ifcn=0
         do i=1,nfxt
            if(abr(i).eq.'XB') then
               ifcn=i
               exit
            endif
         enddo
         CALL DMGFXT(ifcn,IND2)       ! "XB"
         iloc=locd(ind2)
         inum=nwds(ind2)
         CALL DMDLOC('%X_XB',IND2,ISIZ,IPT)
         ilocx=locd(ind2)
         inumx=nwds(ind2)
         if(inum.ne.inumx) then
            call errmsg_exit(
     >         '? ZC X(%XXXB) axis size inconsistency in dmgxot!')
         endif
         datbuf(ilocx:ilocx+inumx-1)=datbuf(iloc:iloc+inum-1)
         DMGLBL(IND2)='%XXXB'
      endif
C
C  %XAZC & %XAZB are usually also copies of X and XB, but they can also
C   be other quantities used in place of X and XB for rplot plotting
C   purposes...
C
C  SEE IF THEY ARE ALREADY THERE
C
      IND1=0
      IND2=0
      CALL DMDLOC('%XAZC',IND1,ISIZ,IPT)
      CALL DMDLOC('%XAZB',IND2,ISIZ,IPT)
C
      IF((IND1.GT.0).AND.(IND2.GT.0)) RETURN
C
C  READ ANY MISSING DATA
C
C  ZONE CTRS
      IF((IND1.EQ.0).AND.(NFX(1).GT.0)) THEN
        CALL DMGFXT(NFX(1),IND1)
        CALL DMGXOT_SPREDM(1,IND1)
        iloc=locd(ind1)
        inum=nwds(ind1)
        call dmdloc('%X_ZC',ind1,ISIZ,IPT)
        ilocx=locd(ind1)
        inumx=nwds(ind1)
        if(inum.ne.inumx) then
           call errmsg_exit('? ZC x axis size inconsistency in dmgxot!')
        endif
        datbuf(ilocx:ilocx+inumx-1)=datbuf(iloc:iloc+inum-1)
        DMGLBL(IND1)='%XAZC'            ! activate
      ENDIF
C  ZONE BDYS
      IF((IND2.EQ.0).AND.(NFX(2).GT.0)) THEN
        CALL DMGFXT(NFX(2),IND2)
        CALL DMGXOT_SPREDM(2,IND2)
        iloc=locd(ind2)
        inum=nwds(ind2)
        call dmdloc('%X_ZB',ind2,ISIZ,IPT)
        ilocx=locd(ind2)
        inumx=nwds(ind2)
        if(inum.ne.inumx) then
           call errmsg_exit('ZB ? x axis size inconsistency in dmgxot!')
        endif
        datbuf(ilocx:ilocx+inumx-1)=datbuf(iloc:iloc+inum-1)
        DMGLBL(IND2)='%XAZB'            ! activate
      ENDIF
      if((ind1.eq.0).and.(ind2.eq.0)) then
         call errmsg_exit(' ?dmgxot -- x axis processing failure!')
      endif
C  IF ZONE CTRS ARE DEFINED BUT NOT BDY'S, INTERPOLATE BDY'S
      IF(IND2.NE.0) GO TO 200
      call dmdloc('%X_ZB',ind2,ISIZ,IPT)
      IXZC=LOCD(IND1)
      ixzb=locd(ind2)
C  LOOP OVER TIMES
      inx=nzonex(1)
      DO 110 IT=1,NTR
         IXZCL=IXZC+(IT-1)*INX
         IXZBL=IXZB+(IT-1)*INX
         CALL XINTZB(DATBUF(IXZCL),DATBUF(IXZBL),INX)
 110  CONTINUE
C  DONE
      DMGLBL(IND2)='%XAZB'              ! activate
      GO TO 1000
C
C  IF ZONE BDYS ARE DEFINED BUT NOT CTRS, INTERPOLATE CTRS
C
 200  CONTINUE
      IF(IND1.NE.0) GO TO 1000
      call dmdloc('%X_ZC',ind1,ISIZ,IPT)
      IXZC=LOCD(IND1)
      IXZB=LOCD(IND2)
C  LOOP OVER TIMES
      inx=nzonex(1)
      DO 210 IT=1,NTR
         IXZCL=IXZC+(IT-1)*INX
         IXZBL=IXZB+(IT-1)*INX
         CALL XINTZC(DATBUF(IXZBL),DATBUF(IXZCL),INX)
 210  CONTINUE
C  DONE
      DMGLBL(IND1)='%XAZC'              ! activate
C--------------------------------------------------------
 1000 continue
cxx      Call dmprin(20)
      RETURN
      END
C--------------------------------------------------------------------------
      subroutine dmgxot_spredm(jtypx,ind)
C
C  do a fudge if a function is monotonic non-decreasing:  make it
C  strictly monotonic increasing
C
C
      use datmgr_mod
      use cplotr_mod
C
C-----------------------------------
C
      inx=nzonex(jtypx)
C
      ia0=locd(ind)
C
      zmini=xpatch(jtypx)
C
      idecr=0
      zdecr=0.0
      zmaxd=0.0
C
      ieq=0
C
      iincr=0
      zincr=0.0
      zmaxi=0.0
C
      do it=1,ntr
         ia=ia0+(it-1)*inx
         do ix=1,inx-1
            ianex=ia+ix
            iaprv=ianex-1
            if((datbuf(iaprv)-zmini).gt.datbuf(ianex)) then
      	 idecr=idecr+1
      	 zdiff=datbuf(iaprv)-datbuf(ianex)
      	 zdecr=zdecr+zdiff
      	 zmaxd=max(zmaxd,zdiff)
            else if(datbuf(ianex).gt.datbuf(iaprv)) then
      	 iincr=iincr+1
      	 zdiff=datbuf(ianex)-datbuf(iaprv)
      	 zincr=zincr+zdiff
      	 zmaxi=max(zmaxi,zdiff)
            else
      	 ieq=ieq+1
            endif
         enddo
      enddo
C
      if(idecr.gt.0) then
         zdecr=zdecr/float(idecr)
         if(iincr.gt.0) then
            zincr=zincr/float(iincr)
         else
            zincr=0.0
         endif
         write(lunzer(0),2222) iincr,zincr,zmaxi,idecr,zdecr,zmaxd,ieq
 2222    format(' %dmgxot_spredm:  non-monotonic X axis:'/
     >'  #increasing steps: ',i7,' avg & max steps: ',2(1x,1pe11.4)/
     >'  #decreasing steps: ',i7,' avg & max steps: ',2(1x,1pe11.4)/
     >'  #zero steps:       ',i7)
         return
      endif
C
      if(ieq.eq.0) return    ! already strict mono. increasing
C
C  OK data is monotonic non-decreasing; make it strict mono. increasing.
C
      do it=1,ntr
         ia=ia0+(it-1)*inx
         iaeq1=0
         iaeq2=0
         do ix=1,inx-1
            ianex=ia+ix
            iaprv=ianex-1
            if(datbuf(iaprv).ge.datbuf(ianex)) then
      	 if(iaeq1.eq.0) iaeq1=iaprv
      	 iaeq2=ianex
            endif
         enddo
C
C  #consecutive equal points (bail out if more than 3)
C
         ineq=iaeq2-iaeq1+1
         if(ineq.gt.1) then
            if(ineq.gt.4) go to 999
            if(ineq.ge.inx) go to 999
C
            if(iaeq1.eq.ia) then
C  at left extremum
      	 zd1=datbuf(iaeq2)-(datbuf(iaeq2+1)-datbuf(iaeq2))
      	 zd2=datbuf(iaeq2)
      	 irhs=1
            else if(iaeq2.eq.ianex) then
C  at right extremum
      	 zd1=datbuf(iaeq1)
      	 zd2=datbuf(iaeq1)+(datbuf(iaeq1)-datbuf(iaeq1-1))
      	 irhs=0
            else
C  in middle
      	 ztes1=datbuf(iaeq1)-datbuf(iaeq1-1)
      	 ztes2=datbuf(iaeq2+1)-datbuf(iaeq2)
      	 if(ztes1.gt.ztes2) then
      	    irhs=1
      	    zd1=datbuf(iaeq1-1)
      	    zd2=datbuf(iaeq1)
      	 else
      	    irhs=0
      	    zd1=datbuf(iaeq2)
      	    zd2=datbuf(iaeq2+1)
      	 endif
            endif
C
C  find & apply patch
C
            zf=5.0e-7
 10         continue
            zf=2.0*zf
            zadd=zf*(zd2-zd1)
            if(irhs.eq.1) then
C  patch from the right
      	 do ial=iaeq2-1,iaeq1
      	    zold=datbuf(ial+1)
      	    znew=zold-zadd
      	    if(znew.lt.zold) then
      	       datbuf(ial)=znew
      	    else
      	       go to 10  ! more spread
      	    endif
      	 enddo
            else
C  patch from the left
      	 do ial=iaeq1+1,iaeq2
      	    zold=datbuf(ial-1)
      	    znew=zold+zadd
      	    if(znew.gt.zold) then
      	       datbuf(ial)=znew
      	    else
      	       go to 10  ! more spread
      	    endif
      	 enddo
            endif
C
         endif  ! ineq
C
      enddo
C
      write(lunzer(0),1001)
 1001 format(' %dmgxot_spredm:  ',
     >    'data patched to be strict monotonic increasing.')
      go to 1000
C
 999  continue
      write(lunzer(0),1002) ineq
 1002 format(' %dmgxot_spredm:  ',
     >    'non monotonic data patch, ineq=',i5)
C
 1000 continue
      return
      end
C
C******************** END FILE DMGXOT.FOR ; GROUP PLDMGR ******************
C
C  remove current %XAZC and %XAZB definitions as needed.
C  rename workspaces to mark them "empty".
C
      subroutine dmgxot_cleanup
 
      use datmgr_mod
      use cplotr_mod
 
      CALL DMDLOC('%XAZC',IND,ISIZ,IPT)
      IF(IND.GT.0) dmglbl(ind)='%X_ZC'
      CALL DMDLOC('%XAZB',IND,ISIZ,IPT)
      IF(IND.GT.0) dmglbl(ind)='%X_ZB'
 
      CALL DMDLOC('%XXXC',IND,ISIZ,IPT)
      IF(IND.GT.0) dmglbl(ind)='%X_XC'
      CALL DMDLOC('%XXXB',IND,ISIZ,IPT)
      IF(IND.GT.0) dmglbl(ind)='%X_XB'
 
      return
      end
