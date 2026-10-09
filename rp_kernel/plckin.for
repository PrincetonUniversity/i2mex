C                                                                   PLCKIN.for
C--------------------------------------------------------------
C  PLCKIN -- CHECK INPUT OPERANDS TO RPLOT USER CALCULATION
C
C    DETERMINE X AXIS OF OUTPUT FUNCTION
C    PROMPT USER IF X AXIS CANNOT BE DETERMINED FROM OPERANDS
C    Called by PLCFXT.
C
 
      SUBROUTINE PLCKIN(ISTAT, Iscal, IER)
 
C	    Last change:
C               8/04/99  DMC  LFNDFUNC moved to RPCALC COMMON
C	        1/28/94  TBT  Commented out InxMr
C	        4/02/93  TbT  Allow major radius and Zone Centered/boundary.
C	       11/23/92  TBT  Allow scalar buffer - Istat=-1.
C                             Added Iscal argument.
C	        1/16/91  TBT  added check of LXOUT before calling PLCECH.
C	       10/05/90  TBT  Added LFNDFUNC to the argument list.
C		9/28/90  TBT  Changed check for X-axis conflict since that
C                             is now done in PLCFNXCT.
C
C	ARGUMENTS:
C	    INPUT:  LFNDFUNC - Logical. True if one of the 13 new functions
C		               was found during the parse. Defined in PLCPARSE.
C		    ISTAT    - I*4 type of X-axis of old accumulator.
C                   Iscal    - I*4 pointer to where in ABT(I) is name $TEMP.
C                              This is where accumulator data is stored when
C                              the accumulator is a scalar.
C	    OUTPUT:
C		    ISTAT    - I*4 type of X axis of the new accumulator.
C                   IER      - error code if an error occurred.
 
C-----------------------------
C
      use datmgr_mod
      use cplotr_mod
      use rpcalc_mod
C
C-----------------------------
C
      CHARACTER*1 ZBUFF(32)
C
      CHARACTER*38 ZLGL
C *LISMAK* DATA SATEMENT MOVED 12-APR-88
C
C 1.  CHECK FOR SELF-REFERENTIAL "$" OPERAND
C 2.  CHECK FOR REFERENCE TO EXISTING RPLOT DATA FUNCTION
C 3.  TRY TO DECODE OPERAND AS A FLOATING POINT NUMBER
C
C
C *LISMAK* 12-APR-88  GENERATED LOGICAL DECLARATIONS
      LOGICAL LXOUT
C *LISMAK* END OF GENERATED LOGICAL DECLARATIONS
C
C => PLCKIN               BEGIN DATA BLOCK
C => LISMAK 12-APR-88  MOVED DATA STATEMENTS
      DATA ZLGL/'ABCDEFGHIJKLMNOPQRSTUVWXYZ$_0123456789'/
C => END DATA BLOCK
C ====================
C
      ILOP=LEN(ZOPND(1))
C
      ILENA=LEN(ABR(1))
C
      IER=0
      ISMAX=0
C
      DO 20 I=1,NOPS
C  NON-BLANK LENGTH OF ARGUMENT
         DO 2 IC=ILOP,1,-1
            IF(ZOPND(I)(IC:IC).NE.' ') THEN
               ILOPI=IC
               GO TO 3
            ENDIF
 2       CONTINUE
         ILOPI=ILOP
 3       CONTINUE
C
         IKOPND(I)=0
         IF((ILOPI.EQ.1).AND.(ZOPND(I)(1:1).EQ.'$')) THEN
            IF(ISTAT.EQ.0) THEN
               CALL ZERMSG('%PLCKIN: "$" NOT YET DEFINED AS OPERAND')
            ELSE
C  THIS OPERAND IS "$" ...
               IKOPND(I)=ISTAT
               If (Istat .eq. -1)  Then
                  Iaopnd(I) = Iscal     ! Accumulator is a scalar.
               Else
                  IAOPND(I)=0           ! Accumulator is a profile.
               End If                   !Istat
            ENDIF
         ELSE
C  CHECK IF FORM IS CORRECT FOR RPLOT FCN NAME
            IF(ILOPI.GT.ILENA) GO TO 10
            ZBUFF(ILENA+1)=' '
            DO 6 IC=1,ILENA
               ZBUFF(IC)=ZOPND(I)(IC:IC)
               IF(ZOPND(I)(IC:IC).EQ.' ') GO TO 6
               IF(INDEX(ZLGL,ZOPND(I)(IC:IC)).EQ.0) GO TO 10
 6          CONTINUE
C  TRY TO DECODE AS RPLOT FCN NAME
            INUM=NFT+NFTX
            IANS=ISCMP0(ZBUFF,ABT,INUM,ILENA,ILEN)
            IF(IANS.NE.0) THEN
               IKOPND(I)=-1
               IAOPND(I)=IANS
            ELSE
               IANS=ISCMP0(ZBUFF,ABR,NFXT,ILENA,ILEN)
               IF(IANS.EQ.0) GO TO 10
C  WE HAVE A PROFILE FUNCTION!  READ INTO MEMORY, GET X AXIS TYPE
C  READ THE DATA AND X AXIS IF APPLICABLE
               CALL DMGFXT(IANS,IND)
               CALL PLPRIO(6,IND)
               IAOPND(I)=IND
               IKOPND(I)=ITYPR(IANS)
            ENDIF
            GO TO 19
C
C  TRY TO DECODE AS A FLOATING POINT NUMBER
 10         CONTINUE
C  no -- plcparse should have detected any numerical constants
CX            ZVOPND(I)=UFDCOD(ZOPND(I),I1U,I2U,ICOD)
CX            IF(ICOD.EQ.0) IKOPND(I)=-2
         ENDIF
C
 19      CONTINUE
         ISMAX=MAX(ISMAX,IKOPND(I))
C
         IF(IKOPND(I).EQ.0) IER=1
C
 20   CONTINUE                          ! END OF OPERAND SCAN LOOP
C
      IF(IER.EQ.1) GO TO 100
C
C  DETERMINE OUTPUT TYPE
C
      IF(ISMAX.LE.0) THEN
         IF (LFNDFUNC) THEN
            ISTAT = 1
         ELSE
            ISTAT = -1
         ENDIF                          ! lfndfunc
      ELSE
 
         If (Ismax .Eq. Ntypmr .And. .Not. Lrmajm) Then
            Call Zermsg('?Plckin - Major radius not defined ')
            Ier=2
            Go to 100
         End If                         ! mjr radius check
 
C     .Check for zone bndry or center consistency will be done in PLCFNXCT.
         IF (ISMAX .Gt. 2)  Then
            DO 50 I=1,NOPS
               IF((IKOPND(I).GT.0).AND.(IKOPND(I).NE.ISMAX)) THEN
C
C            .Allow major radius and zone centered/boundary.
                  If( (.not.NLTRANSP) .Or. (Ismax .Ne. Ntypmr) .Or.
     1               (Ikopnd(I) .Ne. 1 .And. Ikopnd(I) .Ne. 2)) Then
                     CALL ZERMSG('?PLCKIN ')
                     CALL ZERMSG(
     1                  '?PLCKIN - OPERAND X AXIS TYPE CONFLICT')
                     CALL ZERMSG(
     1                  '?PLCKIN   INPUT EQUATION IGNORED, $ UNCHANGED')
                     IER=2
                     GO TO 100
                  End If                ! Ntypmr
               ENDIF
 50         CONTINUE
         End If                         ! Ismax
 
         ISTAT=ISMAX
      ENDIF
C
 100  CONTINUE
C
      RETURN
      END
