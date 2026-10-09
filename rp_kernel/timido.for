C******************** START FILE TIMIDO.FOR ; GROUP PLOTR3 ******************
C-----------------------------------------------------------
C  TIMIDO
C  DO A TIME INTEGRAL
C
      REAL FUNCTION TIMIDO(T,F,N,T1,Z1,T2,Z2)
C
C  T(N) F(N)  FCN BEING INTEGRATED AND TIME VECTOR
C  T1   T2    LIMITS OF TIME INTEGRATION
C  Z1   Z2    INDEX NORM LIMITS T(INT(Z1)) IS GREATEST T(J) .LT. T1
C
C  NOTE CALLER ASSURES (1) T1.GT.T2
C                      (2) T1.GE.T(1)
C                      (3) T2.LE.T(N)
C
      REAL T(N),F(N)
C
      TIMIDO=0.0
C
      I1=IFIX(Z1)+1
      I2=IFIX(Z2)
C
        F1=XINTRP(F,N,Z1,0.0)
        F2=XINTRP(F,N,Z2,0.0)
C
      IF(I1.GT.I2) THEN
C  BOTH TIMES ARE INSIDE THE SAME TIME BIN
        ZDELT=T2-T1
        TIMIDO=0.5*(F1+F2)*ZDELT
C
        RETURN
      ENDIF
C  DEAL WITH TIME BINS IN WHICH THE ENDPTS LIE
        ZI1=0.5*(F(I1)+F1)*(T(I1)-T1)
        ZI2=0.5*(F2+F(I2))*(T2-T(I2))
        TIMIDO=TIMIDO+ZI1+ZI2
C  DEAL WITH INTERIOR TIME BINS (IF ANY)
        ILIM=I2-1
      IF(ILIM.LT.I1) RETURN
C
      DO 10 I=I1,ILIM
        IP=I+1
        TIMIDO=TIMIDO+0.5*(F(I)+F(IP))*(T(IP)-T(I))
 10   CONTINUE
C
      RETURN
      END
C-----------------------------------------------------------
C  TIMIVAR
C  DO A TIME AVERAGED RMS VARIANCE
C
      REAL FUNCTION TIMIVAR(avg,T,F,N,T1,Z1,T2,Z2)
C
C  AVG -- avg value against which variance is computed.
C  T(N) F(N)  FCN BEING INTEGRATED AND TIME VECTOR
C  T1   T2    LIMITS OF TIME INTEGRATION
C  Z1   Z2    INDEX NORM LIMITS T(INT(Z1)) IS GREATEST T(J) .LT. T1
C
C  NOTE CALLER ASSURES (1) T1.GT.T2
C                      (2) T1.GE.T(1)
C                      (3) T2.LE.T(N)
C
      REAL T(N),F(N)
      real :: scal
C
      TIMIVAR=0.0
C
      scal = max(1.0,abs(avg)) ! introduced to avoid risk of floating overflow
C
      I1=IFIX(Z1)+1
      I2=IFIX(Z2)
C
      F1=XINTRP(F,N,Z1,0.0)
      F2=XINTRP(F,N,Z2,0.0)
C
      if(t2.le.t1) then
         timivar=0.0
         return
      endif
C
      IF(I1.GT.I2) THEN
C  BOTH TIMES ARE INSIDE THE SAME TIME BIN
        ZDELT=T2-T1
        TIMIVAR=0.5*( ((F1-avg)/scal)**2 + ((F2-avg)/scal)**2 )*ZDELT
      ELSE
C  DEAL WITH TIME BINS IN WHICH THE ENDPTS LIE
        ZI1=.5*( ((F(I1)-avg)/scal)**2 + ((F1-avg)/scal)**2 )*(T(I1)-T1)
        ZI2=.5*( ((F2-avg)/scal)**2 + ((F(I2)-avg)/scal)**2 )*(T2-T(I2))
        TIMIVAR=TIMIVAR+ZI1+ZI2
C  DEAL WITH INTERIOR TIME BINS (IF ANY)
        ILIM=I2-1
C
        DO I=I1,ILIM
          IP=I+1
          TIMIVAR = TIMIVAR + 0.5*
     >         ( ((F(I)-avg)/scal)**2 + ((F(IP)-avg)/scal)**2 )*
     >         (T(IP)-T(I))
        ENDDO
      ENDIF
C
      TIMIVAR= scal * sqrt(TIMIVAR/(T2-T1))
C
      RETURN
      END
C******************** END FILE TIMIVAR.FOR ; GROUP PLOTR3 ******************
