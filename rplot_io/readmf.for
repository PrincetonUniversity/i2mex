C******************** START FILE READMF.FOR ; GROUP PLOTR4 *************
C---------------------------------------------------------------
C  SUBROUTINE READMF
 
 
      SUBROUTINE READMF(IT,IFCN,F,NF)
 
 
C 	01/28/94 TBT Initialize Zdummy
C    D. MC CUNE 30 JULY 1981
C  READ THE PROFILE VARIATION OF A FCN OUT OF THE MF FILE AT
C  SPECIFIED TIME
C    THIS VERSION SUPPORTS MULTI-RECORD FUNCTIONS.
C
C  INPUTS-- IT = TIME PT. INDEX
C         IFCN = INDEX # OF FUNCTION TO BE READ
C           NF = LENGTH (WORDS) OF FCN TO BE READ
C
C  OUTPUTS- F(NF)  THE PROFILE VARIATION OF THE FCN, READ IN AT THE
C     SPECIFIED TIME.
C
 
C
      use cplotr_mod
      use mfblok_mod
C
      DIMENSION F(NF)
C      DIMENSION ZDUMMY(NR0)

      real,allocatable :: zdummy(:)
C
      logical iblki
      integer :: mfblk
C
C-----------------------------------------------------
      allocate(zdummy(nr0))

      Zdummy(1) = 0.0   ! Initialize tbt
 
      if(lrun_x.eq.0) then
         ITYP=ITYPR(IFCN)
         INF=NZONEX(ITYP)
         INREC=NRECX(ITYP)
         iblki=MFBLKI
      else
         ITYP=ITYPR_X(IFCN,lrun_x)
         INF=NZONEX_X(ITYP,lrun_x)
         INREC=NRECX_X(ITYP,lrun_x)
         iblki=MFBLKI_X(lrun_x)
      endif
C
      IF(INF.LE.NF) GO TO 10
C
      write(lunzer(0),9001)
 9001 FORMAT(' ? PROGRAM ERROR-- INSUFFICIENT ARRAY SPACE PROVIDED'/
     1  '   FOR FCN READ, SUBROUTINE READMF CALLING ARGUMENTS')
      call bad_exit
C
 10   CONTINUE
C-------------------------------------------------------------------
C  READ THE PROFILE AT THE SPECIFIED TIME
C
      IF(iblki) THEN
C------------------------
C  NEW BLOCKED MF FILE
C
C  DATA ADDRESSES, RELATIVE TO BASE BLOCK FOR DATA FCN
        ISTART=(IT-1)*INF+1
        IFIN=IT*INF
C  BLOCK ADDRESSES
        if(lrun_x.eq.0) then
           IBLOK0=MFHDR(IFCN+3)
           iluni=MFLUNI
        else
           IBLOK0=MFHDR_X(IFCN+3,lrun_x)
           iluni=MFLUNI_X(lrun_x)
        endif
        IINC1=(ISTART-1)/IBLKSZ
        IBLOK1=IBLOK0+IINC1
        IINC2=(IFIN-1)/IBLKSZ
        IBLOK2=IBLOK0+IINC2
C  DATA ADDRESSES REL. BLOCK ADDRESSES
        ISTART=ISTART-IBLKSZ*IINC1
        IFIN=IFIN-IBLKSZ*IINC2
C  READ START BLOCK IF NEEDED  ! dmc -- disabled MFBLOK recall for safety
C	  no more ... IF(IBLOK1.NE.MFBLK) THEN
        MFBLK=IBLOK1
        READ(iluni,REC=MFBLK,IOSTAT=IOS) MFDATA
        IF (IOS .NE. 0)  THEN
           write(lunzer(0),9101) iluni, MFBLK, IOS
 9101      FORMAT(/' READMF?: ERROR DURING READ OF UNIT ',
     1	I5/ '   RECORD TRYING TO BE READ = ',
     2	I8/ '   IOSTAT = ', I5/ ' STOPPING')
           call bad_exit
        ENDIF				! IOS
C	  no more ... ENDIF
C  IF TIME VARIATION DOES NOT CROSS BLOCK BDY
        IF(IBLOK1.EQ.IBLOK2) THEN
           DO 12 J=1,INF
      	ID=J+ISTART-1
      	F(J)=MFDATA(ID)
 12        CONTINUE
        ELSE
C  TIME VARIATION CROSSES A BDY
           DO 14 ID=ISTART,IBLKSZ
      	J=ID-ISTART+1
      	F(J)=MFDATA(ID)
 14        CONTINUE
C
C  dmc Jun 1992 -- don't forget to read intermediate blocks!
C
           do iblkl=iblok1+1,iblok2-1
      	read(iluni,rec=iblkl) mfdata
      	do id=1,iblksz
      	   j=j+1
      	   f(j)=mfdata(id)
      	enddo
           enddo
C
           MFBLK=IBLOK2
           READ(iluni,REC=MFBLK) MFDATA
C
           DO 16 ID=1,IFIN
      	J=INF-IFIN+ID
      	F(J)=MFDATA(ID)
 16        CONTINUE
        ENDIF
      ELSE
C------------------------
C  OLD STYLE MF FILE
C
C  CALCULATE FIRST AND LAST RECORD #S OF RECORDS TO READ,
C  AND LEFT OVER DUMMY PTS IN LAST RECORD TO READ IF ANY.
C
         if(lrun_x.eq.0) then
            ifr=nfr
            irofff=nrofff(ifcn)
            izones=nzones
            iluni=MFLUNI
         else
            ifr=nfr_x(lrun_x)
            irofff=nrofff_x(ifcn,lrun_x)
            izones=nzones_x(lrun_x)
            iluni=MFLUNI_X(lrun_x)
         endif
C
         IREC=(IT-1)*(IFR+1)+IROFFF-1
         ILEFT=IZONES*INREC-INF
C  READ LOOP
         IN2=0
C
         DO 30 I=1,INREC
            IN1=IN2+1
            IN2=MIN0(INF,(IN1+IZONES-1))
            IF((I.EQ.INREC).AND.(ILEFT.NE.0)) GO TO 20
C  READ A COMPLETE RECORD INTO FCN BUFFER
            IREC=IREC+1
            READ(iluni,REC=IREC) (F(J),J=IN1,IN2)
            GO TO 30
C  READ A PARTIAL RECORD INTO FCN BUFFER, ROUND OUT WITH DUMMY READS
 20         CONTINUE
            IREC=IREC+1
            READ(iluni,REC=IREC) (F(J),J=IN1,IN2),
     1	 (ZDUMMY(JD),JD=1,ILEFT)
 30      CONTINUE
C
      ENDIF
C-----------------
C DONE
 999  continue
      deallocate(zdummy)
      RETURN
C
      END
C******************** END FILE READMF.FOR ; GROUP PLOTR4 ***************
