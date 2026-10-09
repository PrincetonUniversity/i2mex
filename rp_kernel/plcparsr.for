      subroutine plcparsr(zinput,ier)
C
C  RPLOT calculator parser
C
C  dmc -- extracted from plcfxt.for 4 Aug 1999
C    setup call to PLCPARSE and report errors
C
      use cplotr_mod
      use rpcalc_mod

      character*(*) zinput              ! input line
      integer ier                       ! output completion code  0=OK
C
C--------------------------
C
      ier=0
C
C	  .PUT $ AND ALL ABBREVIATIONS OF SCALARS AND PROFILES INTO A C* ARRAY.
      ALLABB(1) = '$'                   ! $ REPRESENTS ACCUMULATOR, F(X,T)
      INUM = NFT+NFTX                   ! # OF SCALARS (READ IN AND GENERATED)
      ILIM = INUM+NFXT+2
      if(ILIM.gt.MAXABB) then
         call bad_exit(
     >      ' ??PLCPARSR: increase MAXABB parameter in RPCALC COMMON')
      endif
C
      DO IN=1,INUM
         ALLABB(IN+1) = ABT(IN)
      enddo
C
      DO IN=1,NFXT
         ALLABB(IN+INUM+1) = ABR(IN)    ! PROFILE ABBREVIATIONS
      enddo
C
      ALLABB(ILIM) = ';'                ! ; MEANS END OF LIST TO PARSE.
 
 
C	  .CLEAN OUT ARRAYS.
      DO IN=1,MAXRPN
         IRPNLIST(IN) = 0
         INPOSLST(IN) = 0
      enddo
 
      DO IN=1,MAXCON
         ZCONLIST(IN) = 0.D0
      enddo
 
      ilast = len_trim(zinput)
 
      CALL PLCPARSE (ZINPUT(1:ILAST)//';',
                                        ! ; MEANS END OF INPUT LINE.
     1                ALLABB,           ! INPUT ABBREVIATION LIST
     2                IRPNLIST,         ! RETURNED REVERSE POLISH NOTATION LIST
     3                ZCONLIST,         ! RETURNED REAL CONSTANTS LIST.
     4                IER )
 
C
C  the following are now set in PLCPARSE via RPCALC COMMON:
CC     A      INPOSLST,		! RETURNED INPUT LINE POSITION LIST
CC     B      LFNDFUNC,         ! Returned logical T if new fnctns found
C
C	  .ERROR handling.
 
      if(ier.ne.0) then
         call plcerr(ier,zinput,ilast)
         go to 199
      endif
C
C	.FIND WHERE THE TOKENS FOR OPERANDS ARE AND SAVE INFO.
C
      NOPS = 0                          ! # OF OPERANDS FOUND IN EQUATION.
      J = 0                             ! J IS THE CUR. POINTER INTO RPN LIST.
      DO 130 I=1,MAXRPN
         IF (J .GE. MAXRPN)  GO TO 131  ! SEARCH IS DONE (max index reached)
         J = J+1

         if (irpnlist(j) .EQ. 0) go to 131  ! SEARCH IS DONE (end of list)

         IF (IRPNLIST(J) .EQ. 41)  THEN ! 41 SAYS NEXT TOKEN IS OPERAND
            J = J+1
            IF (NOPS .GE. IMAXOP)  THEN ! TOO MANY OPERANDS TBT 9/12/90
               IER = MAX(INPOSLST(J),1) ! MAKE SURE NON-ZERO.
               WRITE(lunzer(0), 8130) IMAXOP
 8130          FORMAT (/' ?PLCFXT: ERROR - MORE THAN',
     1            I4, ' OPERANDS' /
     2            '      INPUT EQUATION IGNORED' /
     3            '      "$" REMAINS THE SAME')
               call plcerr(ier,zinput,ilast)
               go to 199
            END IF                      ! NOPS
 
            NOPS = NOPS + 1
            IOP = IRPNLIST(J)
            IPOSITION(NOPS) = IOP
            ZOPND(NOPS) = ALLABB(IOP)
            !  character position in input string
            IPOS_STR(NOPS) = INPOSLST(J) + 1 - len(trim(zopnd(nops)))
         END IF                         ! TOKEN=41
 
 130  CONTINUE                          ! I LOOP
 131  CONTINUE                          ! 130 loop exit
 
 199  continue
      return
      end
