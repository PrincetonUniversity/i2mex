      subroutine rpctini
C
C  initialize scalar (f(t)) calculator
C
      use datmgr_mod
      use cplotr_mod
      use rpcalc_mod
C
      Character*5 Zabbl
C
      logical fungot
C
C----------------------------------------
C
      Expmax = 37.    ! Machine accuracy for preventing overflow.
      Expmin =-36.    ! Machine accuracy for preventing underflow.
C
      ! Put time into a scalar called TIME
      Call PlfTmk4 ( Time, Time, NTT, 'Time', 'Seconds', 'TIME')
 
C	! Create temporary scalar for holding scalar accumulator.
      If (NFTX .LT. 32) Then
          NFTX = NFTX+1
          Ind = NFT+NFTX
          ABT   (Ind) = '%$TEMP'
          call aordr_add(abt,iordrt,ind)
          LabelT(Ind) = 'Scalar Calculator'
          UnitsT(Ind) = 'None'
          Write (Zabbl, '(''%T'',I3.3)')  Nftx
          Call RP_DMGALO(Ntt, Jloc, 6)
          NWDS(Jloc) = NTT
          DMGLBL(Jloc) = Zabbl
 
          Call DMGFOTX(2, IPT, Ierr)   ! Copy %T001 into $TEMP
 
          NptAcc  = IPT + NTT*(NFT-1) ! Point to $temp in Datbuf
          Nscalar = NFT
      Else
          WRITE(lunzer(0),9901) Nfxt
 9901     FORMAT(' ?rpctini: Error - NFXT too big =', I5 /
     1             '         Scalar calulator init error'/)
          call abortt
      End If ! Nftx
C
C----------------------
C     tbt 6/95
      PltAbb = 'RMAJM'
      Lrmajm = FunGot(IFcnMr)           ! Define major radius type
      Ntypmr = -99
      If (Lrmajm) Ntypmr = Itypr(IFcnMr)
C
      return
      end
