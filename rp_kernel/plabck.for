C-----------------------------------------------------------------------
C  PLABCK -- CHECK THAT INPUT FUNCTION ID DOES NOT ALREADY EXIST
C
      SUBROUTINE PLABCK(ZABR,IER)

      use cplotr_mod
C
      CHARACTER*(*) ZABR
C
      CHARACTER*1 ZBUFF(32)
C
C  ZABR -- INPUT TEST FCN ABBREVIATION
C  IER -- OUTPUT COMPLETION CODE:
C   =0  ID IS UNIQUE
C   =1  ID EXISTS ALREADY OR CONTAINS ILLEGAL CHARACTERS OR IS BLANK
C
      IER=0
C
      LUNT=LUNZER(0)
C
      IF(ZABR.EQ.' ') THEN
        IER=1
        CALL ZERMSG('?PLABCK:  FCN ID IS BLANK')
        GO TO 500
      ENDIF
C
      call idchek0(zabr,ichk,1)
      if(ichk.eq.99) then
         ier=1
         CALL ZERMSG('?PLABCK:  ILLEGAL CHARACTER IN FCN ID: '//zabr)
         GO TO 500
      endif
C
      ILB=LEN(ABT(1))
C
      ila=min(31,len(zabr))
      do i=1,ila
         zbuff(i)=zabr(i:i)
         if(zbuff(i).ne.' ') iblank=0
      enddo
      zbuff(ila+1)=' '
C
      INUM=NFT+NFTX
      IANS=ISCMP0(ZBUFF,ABT,INUM,ILB,ILEN)
      IF(IANS.NE.0) THEN
        IER=1
        WRITE(LUNT,9001) ABT(IANS),LABELT(IANS),UNITST(IANS)
        CALL ZERMSG('?PLABCK:  ID ALREADY IN USE AS SCALAR FUNCTION')
      ENDIF
C
 9001   FORMAT(/' ABBREVIATION:  ',A/' LABEL:  ',A/' UNITS:  ',A)
C
      IANS=ISCMP0(ZBUFF,ABR,NFXT,ILB,ILEN)
      IF(IANS.NE.0) THEN
        IER=1
        WRITE(LUNT,9001) ABR(IANS),LABELR(IANS),UNITSR(IANS)
        CALL ZERMSG('?PLABCK:  ID ALREADY IN USE AS PROFILE FUNCTION')
      ENDIF
C
      IANS=ISCMP0(ZBUFF,ABB,NBAL,ILB,ILEN)
      IF(IANS.NE.0) THEN
        IER=1
        WRITE(LUNT,9001) ABB(IANS),LABELB(IANS),UNITSB(IANS)
        CALL ZERMSG('?PLABCK:  ID ALREADY IN USE AS MULTIGRAPH')
      ENDIF
C
 500  CONTINUE
      RETURN
      END
