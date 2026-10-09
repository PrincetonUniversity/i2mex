C-----------------------------------------------------------------------
C  TDB_EXECSYM -- execute presymmetrization of a profile at a few time pts
C
C  NEW JULY 1991 -- DMC
C   mod dmc May 2005 -- modularization; remove trcom dependency.
C
      SUBROUTINE TDB_EXECSYM(d,
     >     ZDATID,IMAP,ILX,INX,ILF,
     >     ILXSY,INXSY,ILFSY,ILSSY,
     >     IER)
C
      use trdatbuf_obj
      use tdbsub_uts

      IMPLICIT NONE

      type (trdatbuf) :: d

      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ilx,inx,ilf,ilxsy,inxsy,ilfsy,ilssy,ier,imap,it1
      integer irbz1,ilt1,int1
      INTEGER iwork,isiz,j,ioff,it,itm1,itt,ix,ia1,ia2,iout,jm1
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 zx0,zxmax,ztime,zt1,zt2,zfac2,zfac1,zdata
      real*8 :: zrbz,zr1,zr2
!============
      Character*(*) zdatid
C
C  INPUTS
C    d -- the trdatbuf object being operated on
C
C    IMAP=1 --> INPUT PROFILE VS. NORMALIZED MAJOR RADIUS
C    IMAP=2 --> INPUT PROFILE VS. ECE FREQUENCY IN GHZ
C
C    presymmetrization is now done at all (NTIME2) time points; 
C    use of normalized midplane aspect ratio output grid enables this to be
C    done without prior knowledge of flux surface equilibrium geometry.
C
C    ILX -- ADDRESS OF X AXIS OF INPUT PROFILE DATA
C    INX -- NUMBER OF PTS IN X AXIS OF INPUT PROFILE DATA
C    ILF -- ADDRESS OF INPUT PROFILE DATA ITSELF
C
C  OUTPUTS
C    ILXSY -- ADDRESS OF X AXIS OF DATA AFTER SYMMETRIZATION
C    INXSY -- NUMBER OF PTS OF X AXIS OF SYMMETRIZED DATA
C    ILFSY -- ADDRESS OF THE SYMMETRIZED DATA
C    ILSSY -- THE SHIFT IMPLIED BY THE SYMMETRIZATION OF THIS PROFILE
C      THE DATA ITSELF IS PUT INTO THE COMMON ARRAY DATBUF
C
C    IER = 0 --> NORMAL SUCCESSFUL COMPLETION
C    IER = 1 --> ERROR
C
C----------------------------------
C
C  THIS SUBROUTINE CALLS ESSENTIALLY THE SAME CODE THAT HAS LONG BEEN
C  USED IN TRANSP FOR PRESYMMETRIZATION OF DATA (FORMERLY CALLED FROM
C  SUBROUTINE TR2INI FROM AUXVAL).  THE ONION SKIN ALGORITHM IS ENTERED
C  VIA SUBROUTINE PRESYM.  CODE COMMENTS DESCRIBE DETAILS OF THE
C  ALGORITHM; A DESCRIPTION CAN ALSO BE FOUND IN $ TRANSPHELP under
C  $ TRANSPHELP OPER NAMELIST DATA_SYMM ...
C
C  The interface to PRESYM has been modified slightly.  The internal
C  workspaces WORKBUF(LBX(j)..LBX(j)+NBX(j)-1) are used!  These work
C  areas are set up in TR2INI which is called once in subroutine
C  AUXVAL at the start of a TRANSP run.
C
C----------------------------------
      integer :: itest,nzones,nonlin,lunmsg_tdb
C----------------------------------
C
      IER=0
      IF(ILF.LE.0) RETURN      ! DO NOTHING IF THERE IS NO DATA
C
      nonlin = lunmsg_tdb(0)
C
      if(imap.eq.2) then
         if(d%ldatrbz.eq.0) then
            write(nonlin,*) ' ?TDB_EXECSYM -- external magnetic field'
            write(nonlin,*) '  (R*Bz) vs. time must be provided for'
            write(nonlin,*) '  presymmetrization of ECE Temperature.'
            ier=1
            return
         endif
         !  time vector location and size for RBZ
         ilt1=d%ltime1
         int1=d%ntime1
         irbz1=d%ldatrbz
      endif
C
C  PTR TO WORKSPACE (CF CALLER, MAKSYM.FOR)
C
      IWORK=d%LFREE_W
C
C  SIZE OF X AXIS FOR SYMMETRIZED DATA
      INXSY=INX/2+1
C
C  ALLOCATE SPACE FOR X AXIS OF SYMMETRIZED DATA
      CALL TDB_DATALO(d,INXSY,ILXSY)
C
C  ALLOCATE SPACE FOR THE SYMMETRIZED DATA ITSELF
C  ALLOCATE SPACE FOR THE SYMMETRIZATION IMPLIED SHIFT PROFILE
      ISIZ=INXSY*d%NTIME2
      CALL TDB_DATALO(d,ISIZ,ILFSY)
      CALL TDB_DATALO(d,ISIZ,ILSSY)
C
C  DEFINE NEW X VECTOR OF SYMMETRIZED DATA
C  this will be based on minor radius or midplane aspect ratio -- 0 = mag 
C  axis, 1 = plasma boundary, (R2-R1)/(R2+R1) where R1 and R2 are
C  the inner and outer midplane R intercepts, would be aspect ratio (used
C  for ECE data).
C
      ZX0=0.0_R8
      ZXMAX=1.075_R8
C
      DO 20 J=1,INXSY
        IOFF=J-1
        d%DATBUF(ILXSY+IOFF)=ZX0+((ZXMAX-ZX0)*IOFF)/(INXSY-1)
 20   CONTINUE
C
C  NOW SYMMETRIZE THE DATA PROFILE AT EACH TIME IN THE SYMMETRIZATION
C  TIMEBASE
C
      WRITE(NONLIN,1001) ZDATID
 1001 FORMAT(' %TDB_EXECSYM:  PRESYMMETRIZE THE "',A,'" DATA')
C
      DO 60 IT=1,d%NTIME2
C
C  current time
C
        ITM1=IT-1
        ZTIME=d%DATBUF(d%LTIME2+ITM1)
C
C  get midplane intercepts
C
        zr1=d%datbuf(d%lrmp_bdy1+itm1)
        zr2=d%datbuf(d%lrmp_bdy2+itm1)
C
C  COPY THE ORIG DATA X AXIS-- DO FREQUENCY TO RADIUS MAP AND ORDER
C  REVERSAL IF THIS IS ECE DATA
C
        IF(IMAP.EQ.1) THEN
C  COPY
           d%WORKBUF(iwork:iwork+inx-1)=d%DATBUF(ilx:ilx+inx-1)
        ELSE
C  ECE MAP (get vacuum field)
           call tdbsub_lookup(d%datbuf(ilt1:ilt1+int1-1),int1,ztime,
     >          it1,zfac1)
           zrbz=d%datbuf(irbz1+it1-1)*(ONE-zfac1)+
     >          d%datbuf(irbz1+it1)*zfac1
           CALL TDB_SYMXMAP(d,zr1,zr2,zrbz,ILX,INX,IWORK)
        ENDIF
C
C  COPY THE PROFILE DATA AND STORE; REVERSE ORDER IF THIS IS
C  ECE DATA
C
        DO 40 IX=1,INX
           IA1=ILF+(IX-1)*d%NTIME2+ITM1
           ZDATA=d%DATBUF(IA1)
           IF(IMAP.EQ.1) THEN
              IOUT=d%LBX(1)+IX-1
           ELSE
              IOUT=d%LBX(1)+INX-IX
           ENDIF
           d%WORKBUF(IOUT)=ZDATA
 40     CONTINUE
C
C  PRESYMMETRIZE THE DATA
C
        CALL TDB_PRESYM(d,zr1,zr2,IMAP,IWORK,INX,ILXSY,INXSY)
C
C  SAVE THE RESULTS
C
        DO 50 J=1,INXSY
           JM1=J-1
           IOFF=d%NTIME2*JM1+ITM1
C
C  THE DATA ITSELF
           d%DATBUF(ILFSY+IOFF)=d%WORKBUF(d%LBX(7)+JM1)
C
C  THE CORRESPONDING IMPLIED SHIFT PROFILE
           d%DATBUF(ILSSY+IOFF)=d%WORKBUF(d%LBX(8)+JM1)
 50     CONTINUE
C
 59     CONTINUE
 60   CONTINUE
C
      RETURN
      END
