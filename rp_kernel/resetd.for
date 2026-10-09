C******************** START FILE RESETD.FOR ; GROUP PLDMGR ******************
C-------------------------------------
C  RESETD
C
C  COPY 1 NUMBER INTO ARRAY (REALS)
C
      SUBROUTINE RESETD(Y,N,X)
      DIMENSION Y(N)
C
      DO 10 I=1,N
        Y(I)=X
 10   CONTINUE
C
      RETURN
      END
C******************** END FILE RESETD.FOR ; GROUP PLDMGR ******************
