C******************** START FILE TIMIUN.FOR ; GROUP PLOTR3 ******************
C-----------------------------------------------------------
C  TIMIUN
C   TIME INTEGRAL UNITS TRANSFORMS
C
      SUBROUTINE TIMIUN(ZUIN,ZUOUT,IFLAG)
C
C  ZUIN = UNITS OF INTEGRAND
C  ZUOUT = UNITS OF INTEGRAL
C  IFLAG=0 IF TRANSFORM KNOWN, 1 IF NOT KNOWN
C
      PARAMETER (IKNOW=6)
      CHARACTER*(*) ZUIN,ZUOUT
      CHARACTER*10 ZTABL(2,IKNOW)
C
C  TABLE OF KNOWN TRANSFORMS:
C
      DATA ZTABL/'N/CM3/SEC ','N/CM3     ',
     >             'N/CM2/SEC ','N/CM2     ',
     >             'N/SEC     ','N         ',
     >             'WATTS/CM3 ','JLES/CM3  ',
     >             'WATTS/CM2 ','JLES/CM2  ',
     >             'WATTS     ','JOULES    '/
C
      IFLAG=1
      ZUOUT=ZUIN
C
      DO IU=1,IKNOW
C  FIND FIRST AND LAST NONBLANK NON $ CHARACTERS
        CALL TMFLLM(ZUIN,I1I,I2I)
        CALL TMFLLM(ZTABL(1,IU),I1K,I2K)
        IF((I2K-I1K).EQ.(I2I-I1I)) THEN
         IF(ZUIN(I1I:I2I).EQ.ZTABL(1,IU)(I1K:I2K)) THEN
          IFLAG=0
          ZUOUT=ZTABL(2,IU)
          RETURN
         ENDIF
        ENDIF
      ENDDO
C
      RETURN
      END
C******************** END FILE TIMIUN.FOR ; GROUP PLOTR3 ******************
