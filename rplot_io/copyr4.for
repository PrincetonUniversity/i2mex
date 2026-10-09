C******************** START FILE copyr4.FOR ; GROUP PLDMGR ******************
C---------------------------------------------------------------
C  copyr4
C
C  COPY 1 ARRAY INTO ANOTHER (REALS)
C
      SUBROUTINE copyr4(A,B,N)
      DIMENSION A(N),B(N)
C
      DO 10 I=1,N
      B(I)=A(I)
 10   CONTINUE
      RETURN
      END
C******************** END FILE copyr4.FOR ; GROUP PLDMGR ******************
