C******************** START FILE IXCAL0.FOR ; GROUP IXCALC ******************
C----------------------------
C  IXCAL0
C  INTEGRATE ZF1 BY ELEMENTS ZDX, STORE IN ZF2
C  MULTIPLY IN ZSGN; sum of elements stored back in ZF1, destroying
C  original contents.
C
      SUBROUTINE IXCAL0(ZF1,ZF2,ZDX,INX,ZSGN)
C
      REAL ZF1(INX),ZF2(INX),ZDX(INX)
C
      REAL*8 ZSUM,ZSUMD,ZTERM
      REAL*8 :: ZMAX = 1.2345e36
C
      ZSUM=0.0
      ZSUMD=0.0
C
      DO 10 IX=1,INX
         ZSUMD=ZSUMD+ZDX(IX)
         ZTERM=ZF1(IX)*ZDX(IX)*ZSGN
         ZSUM=ZSUM+ZTERM

         IF(ZSUM.gt.zmax) then
            zsum=zmax
         else if(zsum.lt.-zmax) then
            zsum=-zmax
         endif

         ZF2(IX)=ZSUM
         ZF1(IX)=ZSUMD

 10   CONTINUE
      RETURN
      END
C******************** END FILE IXCAL0.FOR ; GROUP IXCALC ******************
