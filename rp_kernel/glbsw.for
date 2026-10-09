C******************** START FILE GLBSW.FOR ; GROUP PLOTR1 ******************
C----------------------------------------------------------------
C  GLBSW
C  SWITCH UNITS LABEL ZLAB IF FOUND IN TABLE ZLSWCH
 
 
      SUBROUTINE GLBSW(ZLSWCH,NSWCH,ZLAB,ISW)
 
 
C	Updated:
C	01/27/94 tbt Added ZZlab=' ' to stop Flint error.
 
      CHARACTER*(*) ZLAB
C
      CHARACTER*10 ZLSWCH(2,NSWCH),ZZLAB,zzswch
C
C---------------------------------
 
      ZZLAB=ZLAB
C
      CALL UUPPER(ZZLAB)
C
      ISW=0
      DO 100 I=1,NSWCH
        IND1=INDEX(ZZLAB,'$')
        IF(IND1.EQ.0) IND1=11
        IND2=INDEX(ZLSWCH(1,I),'$')
        IF(IND2.EQ.0) IND2=11
C  DONT CHECK TRAILING BLANKS
 10     IND1=IND1-1
          IF(IND1.EQ.0) GO TO 20
          IF(ZZLAB(IND1:IND1).EQ.' ') GO TO 10
 20     IND1=MAX0(1,IND1)
C
 30     IND2=IND2-1
          IF(IND2.EQ.0) GO TO 40
          IF(ZLSWCH(1,I)(IND2:IND2).EQ.' ') GO TO 30
 40     IND2=MAX0(1,IND2)
C  CHECK EQUALITY
        IF(IND1.NE.IND2) GO TO 100
          ZZSWCH=ZLSWCH(1,I)
          CALL UUPPER(ZZSWCH)
          IF(ZZLAB(1:IND1).EQ.ZZSWCH(1:IND2)) GO TO 200
C  END OF LLOP
 100  CONTINUE
C
C  FELL THRU:  NO MATCH
      RETURN
C
C  MATCH
C
 200  CONTINUE
      ISW=1
      ZLAB=ZLSWCH(2,I)
      RETURN
      END
C******************** END FILE GLBSW.FOR ; GROUP PLOTR1 ******************
