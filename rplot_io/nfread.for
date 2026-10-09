C******************** START FILE NFREAD.FOR ; GROUP PLOTR1 *************
C==============================================================
C  NFREAD  READ ####NF.PLN FILE (TRANSP GRAPHICS OUTPUT)
C          CONTAINING SCALAR FUNCTIONS OF TIME)
C
C  redone DMC 17 Nov 1997; can be used in context of "primary" or
C  "secondary" run.  See COMMON var. LRUN_X.
C
      SUBROUTINE NFREAD(LUN,IPT,ISIZE,IER)
C
C  IPT-- PTR TO LOCATION AT WHICH TO STORE DATA (don't store if IPT=0)
C  ISIZE-- on input:
C    ISIZE=0 -- need to determine space needed to store data.
C      read data in the order it is found in the file.  This is not
C      time contiguous.
C    ISIZE.gt.0 -- amount of space needed is already known.  This
C      means the number of time points to be read is known; use this
C      info to rearrange data so that it is stored time contiguous
C      and on output, store ISIZE=size of data (#fcns x #times)
C
C  LUN-- LOGICAL UNIT NUMBER ON WHICH TO READ (FILE PRESUMED OPEN)
C
C  IER   COMPLETION CODE  0=SUCCESS
C
C   INCLUDE COMMON BLOCKS ---
      use datmgr_mod
      use cplotr_mod
      use mfblok_mod
C
C------------------------------------------------------------------
C
      real zzbuf(NAXFOT+1)
C
C  START OF EXECUTABLE CODE
C
      IER=0
C
      if(lrun_x.eq.0) then
         if(isize.gt.0) inumt=ntt
         inumf=nft
      else
         if(isize.gt.0) inumt=ntt_x(lrun_x)
         inumf=nft_x(lrun_x)
      endif
C
      if(isize.gt.0) call dmg_texpand(inumt)
C
      ibuf=inumf+1
C
      IREC=0
      IP=IPT-1
C
 10   CONTINUE
C
      IP0=IP+1
      IP=IP+inumf
C
      IREC=IREC+1
C
      READ(LUN,END=100,ERR=1100) (zzbuf(i),i=1,ibuf)
C
      if(irec.ge.ntime) call dmg_texpand(0)

      if(ipt.gt.0) then

         if(lrun_x.eq.0) then
            time(IREC)=zzbuf(1)
         else
            time_x(IREC,lrun_x)=zzbuf(1)
         endif
         if(isize.eq.0) then
            do if=1,inumf
      	 datbuf(IP0+if-1)=zzbuf(1+if)
            enddo
         else
            do if=1,inumf
      	 iadr=ipt+(if-1)*inumt+(irec-1)
      	 datbuf(iadr)=zzbuf(1+if)
            enddo
         endif
      endif
C
      INT0=NTIME
      IF(IREC.LT.INT0) GO TO 10
C
      WRITE(lunzer(0),2001) INT0
 2001 FORMAT(//' % NF FILE NUMBER OF TIME POINTS EXCEEDS PROGRAM'/
     >         '   CAPACITY OF ',I4,' TIME POINTS; EXCESS IGNORED'//)
      IREC=INT0+1
C
 100  CONTINUE
      IREC=IREC-1
      ISIZE=IREC*INUMF
      if(lrun_x.eq.0) then
         ntt=irec
         if(ipt.gt.0) then
            CALL TIMCK1(TIME,NTT)	! DMC 2/89 -- CHECK MONOTONICITY
         endif
      else
         ntt_x(lrun_x)=irec
         if(ipt.gt.0) then
            call timck1(time_x(1,lrun_x),ntt_x(lrun_x))
         endif
      endif
C
C  REWIND FILE-- IN CASE IT NEEDS TO BE REREAD LATER
      REWIND LUN
      RETURN
C
C  ERRORS *****
C
 1100 CONTINUE
      WRITE(lunzer(0),2004)
 2004 FORMAT(' ?NFREAD:  ERROR READING FILE  NF.PLN  DATA')
      IER=-2
C
      RETURN
      	END
C******************** END FILE NFREAD.FOR ; GROUP PLOTR1 ***************
