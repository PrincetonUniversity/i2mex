C******************** START FILE TMFLLM.FOR ; GROUP PLOTR3 ******************
C
      SUBROUTINE TMFLLM(ZU,I1,I2)
      CHARACTER*(*) ZU
C
      ILU=LEN(ZU)
C
      DO 10 I=1,ILU
        IF((ZU(I:I).NE.' ').AND.(ZU(I:I).NE.'$')) GO TO 15
 10   CONTINUE
      I=ILU
 15   CONTINUE
      I1=I
C
      DO 20 I=ILU,I,-1
        IF((ZU(I:I).NE.' ').AND.(ZU(I:I).NE.'$')) GO TO 25
 20   CONTINUE
      I=I1
 25   CONTINUE
      I2=I
C
      RETURN
C
      END
C******************** END FILE TMFLLM.FOR ; GROUP PLOTR3 ******************
