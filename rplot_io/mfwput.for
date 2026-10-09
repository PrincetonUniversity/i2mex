C----------------------------------------------------------------------
C  MFWPUT -- WRITE 512 BYTE BLOCKED MF DATA FILE FOR RPLOT
C
C  DMC APR 1987
C
      SUBROUTINE MFWPUT(LUN,IRECMF,ZDATA,ILEN)
C
      use cplotr_mod
      use mfblok_mod
C
      REAL ZDATA(ILEN)
C
      DO JJ=1,IBLKSZ
        MFDUM(JJ)=0.0
      ENDDO
C
      DO II=1,ILEN,IBLKSZ
        IWD1=II
        IWD2=II+IBLKSZ-1
        IF(IWD2.LE.ILEN) THEN
          IRECMF=IRECMF+1
          WRITE(LUN,REC=IRECMF) (ZDATA(J),J=IWD1,IWD2)
        ELSE
          IREM=IWD2-ILEN
          IRECMF=IRECMF+1
          WRITE(LUN,REC=IRECMF) (ZDATA(J),J=IWD1,ILEN),
     >        (MFDUM(JJ),JJ=1,IREM)
        ENDIF
      ENDDO
C
      RETURN
      END
