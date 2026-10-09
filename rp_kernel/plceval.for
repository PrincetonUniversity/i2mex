      subroutine plceval(istat,iwrk2,zinput,ilast,ier,iwarn,lshift)
C
      use datmgr_mod
      use cplotr_mod
      use rpcalc_mod
C
C  dmc Aug 1999 -- RPLOT calculator evaluation master routine
C    extracted from old plcfxt subroutine.
C
C  this code has no UREAD or TRGRAF dependence and could be called
C  by another application...
C
C  passed input:
      integer istat                     ! type code for expected result
C        -1:  f(t); 1,2,3,... various types of profile evolutions f(x,t)
      integer iwrk2                     ! workspace address
      character*(*) zinput              ! character expression being evaluated
      integer ilast                     ! last non-blank in zinput
C
C  passed output:
      integer ier                       ! return code 0= OK
      integer iwarn                     ! warning code 0= no warning
      logical LSHIFT                    ! flag if zone/bdy shift was imposed
C
C---------------------------------------
C  local...
C
      Integer IkopTemp(ImaxOp)  ! Temporary storage of Ikopnd
C
      INTEGER IARGLIST (NAXFOT+NAXFXT)  ! USED AS A POINTER.
				       ! # OF THE OPERAND (1-NOPS) IN WRKSPC3.
                                       ! VALUE ARGLIST(I) IS FOR ALLABB(I)
                                       ! NOTE: AFTER PARSING - ONLY THOSE
                                       !  VALUES NEEDED WILL BE PUT IN ARGLIST.
CC    INTEGER      IPOSITION(IMAXOP) ! POSITION IN IRPNLIST OF EACH OPERAND.
C  in RPCALC now
C
      DOUBLE PRECISION RESULT(NR0)
C
      LOGICAL ILOCERR                  ! LOCAL ERROR CODE.
      logical idlock
C
C---------------------------------------
C
      iwarn = 0
      ier = 0
C
      IND3   = 0                        ! INITIALIZE
      Ind5   = 0

      call rpcaldat_exec
 
      If (Istat .eq. -1)  Then
         Inx = 1                        ! Scalar: f(t)
         intl = ntt
      Else
         INX = NZONEX(ISTAT)            ! TBT
         intl = ntr
      End If                            ! Istat
 
		! define new work space for holding a row (INX) of data for
		!        each operand (NOPS operands).
      ISIZ3 = INX * NOPS
      CALL RP_DMGALO (ISIZ3, IND3, 10)    ! allocate space
      DMGLBL(IND3) = '%WRK3'            ! label space
      NWDS(IND3)   = ISIZ3              ! set size of workspace
      CALL DMDLOC ('%WRK3', IND3, ISIZ3, IWRK3)
				                     ! find place in buffer
 
 
		! define new work space for holding contents of accumulator
		!        Accumulator may change size (type 2 * type 7)
                !        This prevents overwriting accumulator during calc.
      ISIZ5 = INX * Intl
      CALL RP_DMGALO (ISIZ5, IND5, 10)    ! allocate space
      DMGLBL(IND5) = '%WRK5'            ! label space
      NWDS(IND5)   = ISIZ5              ! set size of workspace
      CALL DMDLOC ('%WRK5', IND5, ISIZ5, IWRK5)
				                     ! find place in buffer
 
      CALL PLC_HAND                     ! ZERO OUT TRAP HANDLER
C		Write(luntrm, 9111) intl, nops, ntypmr, inx, istat
 9111 Format( //' 9111 intl, nops, ntypmr, inx, istat ' /
     1	     5I10 )
 
 
      Do 140 Io=1,Nops
         Ikoptemp(Io) = Ikopnd(Io)      ! Save off to use every time.
 140  Continue
 
C		Write(luntrm, 9112) (ikopnd(io), io=1,nops)
 9112 Format( ' ikopnd = ', ((5I10) /))
 
      DO 180 IT=1,intl
C		    Write(luntrm, 9113) it
 9113    Format(' New it loop ', i5)
 
         Do 145  IO=1,Nops
            Ikopnd(Io) = Ikoptemp(IO)   ! Restore original each time
 145     Continue                       !
 
C		    Write(luntrm, 9111) intl, nops, ntypmr, inx, istat
C		    Write(luntrm, 9112) (ikopnd(io), io=1,nops)
 
         DO 160  IO=1,NOPS
 
C		        Write(luntrm, 9114) io
 9114       Format( ' New io Loop     ', I5)
 
            DO 150 IX =1,INX
C		  	    IADL = (IWRK2-1) + (IT-1)*INX + IX  tbt
C                                     do above calculation in PLCCOP
 
			    ! ONLY LOAD NEEDED VALUES.
 
               IADL3 = (IWRK3-1) + (IO-1)*INX + IX
               IPOS = IPOSITION(IO)
               IARGLIST(IPOS) = IO
               DATBUF(IADL3) =PLCCOP(Istat, Io,
     1              IT,IX,INX,Iwrk2,
     2              Ntypmr, nptacc)
 150        Continue
 
            If((Istat .eq. Ntypmr .And. Ikopnd(IO) .Eq. 1) .Or.
     1         (Istat .eq. Ntypmr .And. Ikopnd(IO) .Eq. 2)) Then
 
               Ikopnd(IO) = Ntypmr
 
C                        PLCCOP explodes type 1 & 2 profiles
C                                (Zone Centered & Boundary)
C                        to major radius (Ntypmr) if any profile is type Ntypmr
            End If                      ! Istat=Ntypmr
 
 160     CONTINUE                       ! END DO IO
 
		    ! Calculate RESULT
         CALL PLCFNXCT (RESULT,   IRPNLIST,
     A                             ZCONLIST, IARGLIST, PREC,
     1                             DATBUF(IWRK3), INX, IT,
     2                             ISTAT,  ZINPUT,
     3				   ILOCERR, IER, LSHIFT)
 
         IF (ILOCERR) THEN
            call plcerr(ier,zinput,ilast)
            go to 1999
         END IF                         ! ilocerr
 
         DO 170 IX=1,INX
            IADL = (IWRK5-1) + (IT-1)*INX + IX
            DATBUF(IADL) = RESULT(IX)
 
C			    If (it .eq. 1 .or. it .eq. 5) Then
C                               Write(luntrm,9212) datbuf(iadl), iadl,
C     1                               io,it,ix,inx,istat,ntypmr,iwrk2
C 9212			       Format( ' Datbuf: ', E15.5, 9I6)
C                            End If
 170     CONTINUE                       ! IX LOOP
 180  CONTINUE                          ! IT LOOP
 
C               Copy temporary area into Accumulator
      Do 1800 It=1,Intl
         DO 1700 IX=1,INX
            IADL2 = (IWRK2-1) + (IT-1)*INX + IX
            IADL5 = (IWRK5-1) + (IT-1)*INX + IX
            DATBUF(IADL2) = Datbuf(Iadl5)
 1700    CONTINUE                       ! IX LOOP
 1800 CONTINUE                          ! IT LOOP
	
C----------------------- TBT
 
C	.Calculation OK.
      PFLAB = Zinput      ! Set label on bottom = Last expression.
      PFLABU= 'Unknown'
C
C	.PRINT OUT ERROR MESSAGES AND ZERO OUT RESULTS WITH ILLEGAL VALUE.
      IPTS = NTR*INX
      CALL PLC_SMRY(DATBUF(IWRK2),IPTS,ZINPUT(1:ILAST),IWARN)
 
 1999 continue
 
      idlock = no_delete
      no_delete = .FALSE.
 
      IF (IND3 .NE. 0)  CALL DMIDEL (IND3)	 ! delete workspace 3
      IND3 = 0
 
      IF (IND5 .NE. 0)  CALL DMIDEL (IND5) ! delete workspace 5
      IND5 = 0
 
      no_delete = idlock
C
C  NOW RESET PRIORITY OF ITEMS STORED FOR THIS CALCLATION TO NORMAL...
C
      DO 500 I=1,NOPS
         IF(IKOPND(I).GT.0) THEN
            IND=IAOPND(I)
            CALL PLPRIO(5,IND)
         ENDIF
 500  CONTINUE
C
      return
      end
