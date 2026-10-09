!
!  "R" routines from fftpack, gathered together, names prepended with "R4"
!      REAL precision -- dmc Aug 2009; implicit none imposed
!
!      IMSL fft replacement effort
!
!      test routine included....
!
!--------------------------------------------
      subroutine R4RFFTI (N,WSAVE)
        use iso_c_binding, only: fp => c_double
        implicit none
        integer ::    N
        real :: WSAVE(*)

        if (N .EQ. 1) return
        CALL R4RFFTI1 (N,WSAVE(N+1),WSAVE(2*N+1))
        return
      end subroutine R4RFFTI

!--------------------------------------------
      subroutine R4RFFTF (N,R,WSAVE)
        use iso_c_binding, only: fp => c_double
        implicit none

        integer ::    N
        real  :: R(*),WSAVE(*)

        if (N .EQ. 1) return
        CALL R4RFFTF1 (N,R,WSAVE,WSAVE(N+1),WSAVE(2*N+1))
        
        return
      end subroutine R4RFFTF
      
!--------------------------------------------
      subroutine R4RFFTB (N,R,WSAVE)
        use iso_c_binding, only: fp => c_double
        implicit none

        integer ::    N
        real :: R(*),WSAVE(*)
        
        if (N .EQ. 1) return
        CALL R4RFFTB1 (N,R,WSAVE,WSAVE(N+1),WSAVE(2*N+1))
        
        return
      end subroutine R4RFFTB

!--------------------------------------------
      subroutine R4RFFTI1 (N,WA,WifAC)
        use iso_c_binding, only: fp => c_double
        implicit none

        integer :: N
        real :: WA(*),WifAC(*)
!      integer ::    ifAC(*)

!-----------------
        integer ::  nl,nf,i,j,ntry,nq,nr,ib,is,l1,nfm1,ip,ld,l2
        integer :: ido,ipm,k1,ii
        real :: tpi,argh,arg,fi,argld

        integer ::    NTRYH(4)
        DATA NTRYH /4,2,3,5/
!-----------------

        NL = N
        NF = 0
        J = 0
101     J = J+1
        if (J-4.le.0) then
           NTRY = NTRYH(J)
        else
           NTRY = NTRY+2
        end if
104     NQ = NL/NTRY
        NR = NL-NTRY*NQ
        if (NR.ne.0) goto 101
        NF = NF+1
        WifAC(NF+2) = NTRY
        NL = NQ
        if (NTRY .NE. 2) goto 107
        if (NF .EQ. 1) goto 107
        do  I=2,NF
           IB = NF-I+2
           WifAC(IB+2) = WIFAC(IB+1)
        end do
        WifAC(3) = 2
107     if (NL .NE. 1) goto 104
        WifAC(1) = N
        WifAC(2) = NF
        TPI = 6.28318530717959
        ARGH = TPI/FLOAT(N)
        IS = 0
        NFM1 = NF-1
        L1 = 1
      if (NFM1 .EQ. 0) return
      do K1=1,NFM1
         IP = int(WifAC(K1+2))
         LD = 0
         L2 = L1*IP
         Ido = N/L2
         IPM = IP-1
         do J=1,IPM
            LD = LD+L1
            I = IS
            ARGLD = FLOAT(LD)*ARGH
            FI = 0.
            do II=3,IDO,2
               I = I+2
               FI = FI+1.
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
            end do
            IS = IS+Ido
         end do
         L1 = L2
      end do

      return
      end

!--------------------------------------------
      subroutine R4RFFTF1 (N,C,CH,WA,WifAC)
use iso_c_binding, only: fp => c_double
      implicit none

      integer :: N
      real  ::     CH(*)      ,C(*)       ,WA(*),  WifAC(*)
!      integer :: ifAC(*)

!-----------------
      integer :: nf,na,l2,iw,k1,kh,ip,l1,ido,idl1,ix2,ix3,ix4,i
!-----------------

      NF = INT(WifAC(2))
      NA = 1
      L2 = N
      IW = N
      do K1=1,NF
         KH = NF-K1
         IP = int(WifAC(KH+3))
         L1 = L2/IP
         Ido = N/L2
         IDL1 = Ido*L1
         IW = IW-(IP-1)*Ido
         NA = 1-NA
         if (IP .NE. 4) goto 102
         IX2 = IW+Ido
         IX3 = IX2+Ido
         if (NA .NE. 0) goto 101
         CALL R4RADF4 (Ido,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         goto 110
  101    CALL R4RADF4 (Ido,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
         goto 110
  102    if (IP .NE. 2) goto 104
         if (NA .NE. 0) goto 103
         CALL R4RADF2 (Ido,L1,C,CH,WA(IW))
         goto 110
  103    CALL R4RADF2 (Ido,L1,CH,C,WA(IW))
         goto 110
  104    if (IP .NE. 3) goto 106
         IX2 = IW+Ido
         if (NA .NE. 0) goto 105
         CALL R4RADF3 (Ido,L1,C,CH,WA(IW),WA(IX2))
         goto 110
  105    CALL R4RADF3 (Ido,L1,CH,C,WA(IW),WA(IX2))
         goto 110
  106    if (IP .NE. 5) goto 108
         IX2 = IW+Ido
         IX3 = IX2+Ido
         IX4 = IX3+Ido
         if (NA .NE. 0) goto 107
         CALL R4RADF5 (Ido,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         goto 110
  107    CALL R4RADF5 (Ido,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         goto 110
  108    if (Ido .EQ. 1) NA = 1-NA
         if (NA .NE. 0) goto 109
         CALL R4RADFG (Ido,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         NA = 1
         goto 110
  109    CALL R4RADFG (Ido,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
         NA = 0
  110    L2 = L1
      end do
       if (NA .EQ. 1) return
      do I=1,N
         C(I) = CH(I)
      end do

      return
      end

!--------------------------------------------
      subroutine R4RFFTB1 (N,C,CH,WA,WifAC)
use iso_c_binding, only: fp => c_double
      implicit none

      integer :: N
      real ::      CH(*)      ,C(*)       ,WA(*),WifAC(*)
!      integer :: ifAC(*)

!-----------------
      integer :: nf,na,l2,iw,k1,kh,ip,l1,ido,idl1,ix2,ix3,ix4,i
!-----------------

      NF =int(WifAC(2))
      NA = 0
      L1 = 1
      IW = 1
      do K1=1,NF
         IP = int(WifAC(K1+2))
         L2 = IP*L1
         Ido = N/L2
         IDL1 = Ido*L1
         if (IP .NE. 4) goto 103
         IX2 = IW+Ido
         IX3 = IX2+Ido
         if (NA .NE. 0) goto 101
         CALL R4RADB4 (Ido,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         goto 102
  101    CALL R4RADB4 (Ido,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         goto 115
  103    if (IP .NE. 2) goto 106
         if (NA .NE. 0) goto 104
         CALL R4RADB2 (Ido,L1,C,CH,WA(IW))
         goto 105
  104    CALL R4RADB2 (Ido,L1,CH,C,WA(IW))
  105    NA = 1-NA
         goto 115
  106    if (IP .NE. 3) goto 109
         IX2 = IW+Ido
         if (NA .NE. 0) goto 107
         CALL R4RADB3 (Ido,L1,C,CH,WA(IW),WA(IX2))
         goto 108
  107    CALL R4RADB3 (Ido,L1,CH,C,WA(IW),WA(IX2))
  108    NA = 1-NA
         goto 115
  109    if (IP .NE. 5) goto 112
         IX2 = IW+Ido
         IX3 = IX2+Ido
         IX4 = IX3+Ido
         if (NA .NE. 0) goto 110
         CALL R4RADB5 (Ido,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         goto 111
  110    CALL R4RADB5 (Ido,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         goto 115
  112    if (NA .NE. 0) goto 113
         CALL R4RADBG (Ido,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         goto 114
  113    CALL R4RADBG (Ido,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114    if (Ido .EQ. 1) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*Ido
      end do
      if (NA .EQ. 0) return
      do I=1,N
         C(I) = CH(I)
      end do

      return
    end subroutine R4RFFTB1

!--------------------------------------------
    subroutine R4RADF2 (Ido,L1,CC,CH,WA1)
      use iso_c_binding, only: fp => c_double
      implicit none

      integer :: Ido,L1
      REAL       CH(Ido,2,L1)           ,CC(IDO,L1,2)           , &
                      WA1(*)

!-----------------
      integer :: k,ic,i,idp2
      real :: tr2,ti2
!-----------------

      do K=1,L1
         CH(1,1,K) = CC(1,K,1)+CC(1,K,2)
         CH(Ido,2,K) = CC(1,K,1)-CC(1,K,2)
      end do
      if (Ido-2.lt.0) return
      if (Ido-2.gt.0) then      !107,105,102
         IDP2 = Ido+2
         do K=1,L1
            do I=3,IDO,2
               IC = IDP2-I
               TR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
               TI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
               CH(I,1,K) = CC(I,K,1)+TI2
               CH(IC,2,K) = TI2-CC(I,K,1)
               CH(I-1,1,K) = CC(I-1,K,1)+TR2
               CH(IC-1,2,K) = CC(I-1,K,1)-TR2
            end do
         end do
         if (MOD(Ido,2) .EQ. 1) return
      end if
      do K=1,L1
         CH(1,2,K) = -CC(Ido,K,2)
         CH(Ido,1,K) = CC(IDO,K,1)
      end do
      return
      end

!--------------------------------------------
      subroutine R4RADF3 (Ido,L1,CC,CH,WA1,WA2)
use iso_c_binding, only: fp => c_double
      implicit none

      integer :: Ido,L1
      REAL       CH(Ido,3,L1)           ,CC(IDO,L1,3)           , &
                      WA1(*)     ,WA2(*)

!-----------------
      integer :: i,ic,k,idp2
      real :: cr2,ci2,dr2,di2,dr3,di3,tr2,ti2,tr3,ti3
      real :: TAUR,TAUI
!-----------------
      DATA TAUR,TAUI /-.5,.866025403784439/

      do K=1,L1
         CR2 = CC(1,K,2)+CC(1,K,3)
         CH(1,1,K) = CC(1,K,1)+CR2
         CH(1,3,K) = TAUI*(CC(1,K,3)-CC(1,K,2))
         CH(Ido,2,K) = CC(1,K,1)+TAUR*CR2
      end do
      if (Ido .EQ. 1) return
      IDP2 = Ido+2
      do K=1,L1
         do I=3,IDO,2
            IC = IDP2-I
            DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
            DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
            CR2 = DR2+DR3
            CI2 = DI2+DI3
            CH(I-1,1,K) = CC(I-1,K,1)+CR2
            CH(I,1,K) = CC(I,K,1)+CI2
            TR2 = CC(I-1,K,1)+TAUR*CR2
            TI2 = CC(I,K,1)+TAUR*CI2
            TR3 = TAUI*(DI2-DI3)
            TI3 = TAUI*(DR3-DR2)
            CH(I-1,3,K) = TR2+TR3
            CH(IC-1,2,K) = TR2-TR3
            CH(I,3,K) = TI2+TI3
            CH(IC,2,K) = TI3-TI2
         end do
      end do

      return
      end

!--------------------------------------------
      subroutine R4RADF4 (Ido,L1,CC,CH,WA1,WA2,WA3)
use iso_c_binding, only: fp => c_double
      implicit none

      integer :: Ido,L1
      REAL       CC(Ido,L1,4)           ,CH(IDO,4,L1)           , &
                      WA1(*)     ,WA2(*)     ,WA3(*)

!-----------------
      integer :: k,idp2,i,ic
      real :: tr1,ti1,tr2,ti2,cr2,ci2,cr3,ci3,cr4,ci4
      real :: tr3,ti3,tr4,ti4
      real :: hsqt2
!-----------------
      DATA HSQT2 /.7071067811865475/

      do K=1,L1
         TR1 = CC(1,K,2)+CC(1,K,4)
         TR2 = CC(1,K,1)+CC(1,K,3)
         CH(1,1,K) = TR1+TR2
         CH(Ido,4,K) = TR2-TR1
         CH(Ido,2,K) = CC(1,K,1)-CC(1,K,3)
         CH(1,3,K) = CC(1,K,4)-CC(1,K,2)
      end do
      if (Ido-2.lt.0) return
      if (Ido-2.gt.0) then
         IDP2 = Ido+2
         do K=1,L1
            do I=3,IDO,2
               IC = IDP2-I
               CR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
               CI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
               CR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
               CI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
               CR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)
               CI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)
               TR1 = CR2+CR4
               TR4 = CR4-CR2
               TI1 = CI2+CI4
               TI4 = CI2-CI4
               TI2 = CC(I,K,1)+CI3
               TI3 = CC(I,K,1)-CI3
               TR2 = CC(I-1,K,1)+CR3
               TR3 = CC(I-1,K,1)-CR3
               CH(I-1,1,K) = TR1+TR2
               CH(IC-1,4,K) = TR2-TR1
               CH(I,1,K) = TI1+TI2
               CH(IC,4,K) = TI1-TI2
               CH(I-1,3,K) = TI4+TR3
               CH(IC-1,2,K) = TR3-TI4
               CH(I,3,K) = TR4+TI3
               CH(IC,2,K) = TR4-TI3
            end do
         end do
         if (MOD(Ido,2) .EQ. 1) return
      end if
      do K=1,L1
         TI1 = -HSQT2*(CC(Ido,K,2)+CC(IDO,K,4))
         TR1 = HSQT2*(CC(Ido,K,2)-CC(IDO,K,4))
         CH(Ido,1,K) = TR1+CC(IDO,K,1)
         CH(Ido,3,K) = CC(IDO,K,1)-TR1
         CH(1,2,K) = TI1-CC(Ido,K,3)
         CH(1,4,K) = TI1+CC(Ido,K,3)
      end do
      return
      end

!--------------------------------------------
      subroutine R4RADF5 (Ido,L1,CC,CH,WA1,WA2,WA3,WA4)
use iso_c_binding, only: fp => c_double
      implicit none

      integer :: Ido,L1
      REAL       CC(Ido,L1,5)           ,CH(IDO,5,L1)           , &
                      WA1(*)     ,WA2(*)     ,WA3(*)     ,WA4(*)

!-----------------
      integer :: i,k,idp2,ic
      real :: tr11,ti11,tr12,ti12
      real :: cr2,cr3,cr4,cr5,ci2,ci3,ci4,ci5
      real :: dr2,dr3,dr4,dr5,di2,di3,di4,di5
      real :: tr2,tr3,tr4,tr5,ti2,ti3,ti4,ti5
!-----------------
      DATA TR11,TI11,TR12,TI12 /.309016994374947,.951056516295154, &
      -.809016994374947,.587785252292473/

      do K=1,L1
         CR2 = CC(1,K,5)+CC(1,K,2)
         CI5 = CC(1,K,5)-CC(1,K,2)
         CR3 = CC(1,K,4)+CC(1,K,3)
         CI4 = CC(1,K,4)-CC(1,K,3)
         CH(1,1,K) = CC(1,K,1)+CR2+CR3
         CH(Ido,2,K) = CC(1,K,1)+TR11*CR2+TR12*CR3
         CH(1,3,K) = TI11*CI5+TI12*CI4
         CH(Ido,4,K) = CC(1,K,1)+TR12*CR2+TR11*CR3
         CH(1,5,K) = TI12*CI5-TI11*CI4
      end do
      if (Ido .EQ. 1) return
      IDP2 = Ido+2
      do K=1,L1
         do I=3,IDO,2
            IC = IDP2-I
            DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
            DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
            DR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)
            DI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)
            DR5 = WA4(I-2)*CC(I-1,K,5)+WA4(I-1)*CC(I,K,5)
            DI5 = WA4(I-2)*CC(I,K,5)-WA4(I-1)*CC(I-1,K,5)
            CR2 = DR2+DR5
            CI5 = DR5-DR2
            CR5 = DI2-DI5
            CI2 = DI2+DI5
            CR3 = DR3+DR4
            CI4 = DR4-DR3
            CR4 = DI3-DI4
            CI3 = DI3+DI4
            CH(I-1,1,K) = CC(I-1,K,1)+CR2+CR3
            CH(I,1,K) = CC(I,K,1)+CI2+CI3
            TR2 = CC(I-1,K,1)+TR11*CR2+TR12*CR3
            TI2 = CC(I,K,1)+TR11*CI2+TR12*CI3
            TR3 = CC(I-1,K,1)+TR12*CR2+TR11*CR3
            TI3 = CC(I,K,1)+TR12*CI2+TR11*CI3
            TR5 = TI11*CR5+TI12*CR4
            TI5 = TI11*CI5+TI12*CI4
            TR4 = TI12*CR5-TI11*CR4
            TI4 = TI12*CI5-TI11*CI4
            CH(I-1,3,K) = TR2+TR5
            CH(IC-1,2,K) = TR2-TR5
            CH(I,3,K) = TI2+TI5
            CH(IC,2,K) = TI5-TI2
            CH(I-1,5,K) = TR3+TR4
            CH(IC-1,4,K) = TR3-TR4
            CH(I,5,K) = TI3+TI4
            CH(IC,4,K) = TI4-TI3
         end do
      end do

      return
      end

!--------------------------------------------
      subroutine R4RADFG (Ido,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
use iso_c_binding, only: fp => c_double
      implicit none

      integer :: Ido,L1,IP,IDL1

      REAL       CH(Ido,L1,IP)          ,CC(IDO,IP,L1)          , &
                      C1(Ido,L1,IP)          ,C2(IDL1,IP), &
                      CH2(IDL1,IP)           ,WA(*)

!-----------------
      real :: tpi,arg,dcp,dsp,ar1,ai1,ar1h,dc2,ds2,ar2,ai2,ar2h
      integer :: ipph,ipp2,idp2,nbd,ik,i,j,k,is,idij,jc,l,lc,j2,ic
!-----------------
      DATA TPI/6.28318530717959/

      ARG = TPI/FLOAT(IP)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IPPH = (IP+1)/2
      IPP2 = IP+2
      IDP2 = Ido+2
      NBD = (Ido-1)/2
      if (Ido .EQ. 1) goto 119
      do IK=1,IDL1
         CH2(IK,1) = C2(IK,1)
      end do
      do J=2,IP
         do K=1,L1
            CH(1,K,J) = C1(1,K,J)
         end do
      end do
      if (NBD .GT. L1) goto 107
      IS = -Ido
      do J=2,IP
         IS = IS+Ido
         IDIJ = IS
         do I=3,IDO,2
            IDIJ = IDIJ+2
            do K=1,L1
               CH(I-1,K,J) = WA(IDIJ-1)*C1(I-1,K,J)+WA(IDIJ)*C1(I,K,J)
               CH(I,K,J) = WA(IDIJ-1)*C1(I,K,J)-WA(IDIJ)*C1(I-1,K,J)
            end do
         end do
      end do
      goto 111
  107 IS = -Ido
      do J=2,IP
         IS = IS+Ido
         do K=1,L1
            IDIJ = IS
            do I=3,IDO,2
               IDIJ = IDIJ+2
               CH(I-1,K,J) = WA(IDIJ-1)*C1(I-1,K,J)+WA(IDIJ)*C1(I,K,J)
               CH(I,K,J) = WA(IDIJ-1)*C1(I,K,J)-WA(IDIJ)*C1(I-1,K,J)
            end do
         end do
      end do
 111  if (NBD .LT. L1) goto 115
      do 114 J=2,IPPH
         JC = IPP2-J
         do 113 K=1,L1
            do 112 I=3,IDO,2
               C1(I-1,K,J) = CH(I-1,K,J)+CH(I-1,K,JC)
               C1(I-1,K,JC) = CH(I,K,J)-CH(I,K,JC)
               C1(I,K,J) = CH(I,K,J)+CH(I,K,JC)
               C1(I,K,JC) = CH(I-1,K,JC)-CH(I-1,K,J)
  112       continue
  113    continue
  114 continue
      goto 121
  115 do J=2,IPPH
         JC = IPP2-J
         do I=3,IDO,2
            do K=1,L1
               C1(I-1,K,J) = CH(I-1,K,J)+CH(I-1,K,JC)
               C1(I-1,K,JC) = CH(I,K,J)-CH(I,K,JC)
               C1(I,K,J) = CH(I,K,J)+CH(I,K,JC)
               C1(I,K,JC) = CH(I-1,K,JC)-CH(I-1,K,J)
            end do
         end do
      end do
      goto 121
 119  do IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
      end do
 121  do J=2,IPPH
         JC = IPP2-J
         do K=1,L1
            C1(1,K,J) = CH(1,K,J)+CH(1,K,JC)
            C1(1,K,JC) = CH(1,K,JC)-CH(1,K,J)
         end do
      end do
!
      AR1 = 1.
      AI1 = 0.
      do L=2,IPPH
         LC = IPP2-L
         AR1H = DCP*AR1-DSP*AI1
         AI1 = DCP*AI1+DSP*AR1
         AR1 = AR1H
         do IK=1,IDL1
            CH2(IK,L) = C2(IK,1)+AR1*C2(IK,2)
            CH2(IK,LC) = AI1*C2(IK,IP)
         end do
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         do J=3,IPPH
            JC = IPP2-J
            AR2H = DC2*AR2-DS2*AI2
            AI2 = DC2*AI2+DS2*AR2
            AR2 = AR2H
            do IK=1,IDL1
               CH2(IK,L) = CH2(IK,L)+AR2*C2(IK,J)
               CH2(IK,LC) = CH2(IK,LC)+AI2*C2(IK,JC)
            end do
         end do
      end do
      do J=2,IPPH
         do IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+C2(IK,J)
         end do
      end do
!
      if (Ido .LT. L1) goto 132
      do K=1,L1
         do I=1,IDO
            CC(I,1,K) = CH(I,K,1)
         end do
      end do
      goto 135
  132 do I=1,IDO
         do K=1,L1
            CC(I,1,K) = CH(I,K,1)
         end do
      end do
  135 do J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         do K=1,L1
            CC(Ido,J2-2,K) = CH(1,K,J)
            CC(1,J2-1,K) = CH(1,K,JC)
         end do
      end do
      if (Ido .EQ. 1) return
      if (NBD .LT. L1) goto 141
      do J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         do K=1,L1
            do I=3,IDO,2
               IC = IDP2-I
               CC(I-1,J2-1,K) = CH(I-1,K,J)+CH(I-1,K,JC)
               CC(IC-1,J2-2,K) = CH(I-1,K,J)-CH(I-1,K,JC)
               CC(I,J2-1,K) = CH(I,K,J)+CH(I,K,JC)
               CC(IC,J2-2,K) = CH(I,K,JC)-CH(I,K,J)
            end do
         end do
      end do
      return
  141 do J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         do I=3,IDO,2
            IC = IDP2-I
            do K=1,L1
               CC(I-1,J2-1,K) = CH(I-1,K,J)+CH(I-1,K,JC)
               CC(IC-1,J2-2,K) = CH(I-1,K,J)-CH(I-1,K,JC)
               CC(I,J2-1,K) = CH(I,K,J)+CH(I,K,JC)
               CC(IC,J2-2,K) = CH(I,K,JC)-CH(I,K,J)
            end do
         end do
      end do
      return
      end

!--------------------------------------------
      subroutine R4RADB2 (Ido,L1,CC,CH,WA1)
use iso_c_binding, only: fp => c_double
      implicit none

      integer :: Ido,L1
      REAL       CC(Ido,2,L1)           ,CH(IDO,L1,2)           , &
                      WA1(*)

!-----------------
      integer :: i,k,ic,idp2
      real :: tr2,ti2
!-----------------

      do K=1,L1
         CH(1,K,1) = CC(1,1,K)+CC(Ido,2,K)
         CH(1,K,2) = CC(1,1,K)-CC(Ido,2,K)
      end do
      if (Ido-2.lt.0) return
      if (Ido-2.gt.0) then
         IDP2 = Ido+2
         do K=1,L1
            do I=3,IDO,2
               IC = IDP2-I
               CH(I-1,K,1) = CC(I-1,1,K)+CC(IC-1,2,K)
               TR2 = CC(I-1,1,K)-CC(IC-1,2,K)
               CH(I,K,1) = CC(I,1,K)-CC(IC,2,K)
               TI2 = CC(I,1,K)+CC(IC,2,K)
               CH(I-1,K,2) = WA1(I-2)*TR2-WA1(I-1)*TI2
               CH(I,K,2) = WA1(I-2)*TI2+WA1(I-1)*TR2
            end do
         end do
         if (MOD(Ido,2) .EQ. 1) return
      end if
      do K=1,L1
         CH(Ido,K,1) = CC(IDO,1,K)+CC(IDO,1,K)
         CH(Ido,K,2) = -(CC(1,2,K)+CC(1,2,K))
      end do
      return
      end

!--------------------------------------------
      subroutine R4RADB3 (Ido,L1,CC,CH,WA1,WA2)
use iso_c_binding, only: fp => c_double
      implicit none

      integer :: Ido,L1
      REAL       CC(Ido,3,L1)           ,CH(IDO,L1,3)           , &
                      WA1(*)     ,WA2(*)

!-----------------
      integer :: i,k,ic,idp2
      real :: tr2,ti2,cr2,ci2,cr3,ci3,dr2,di2,dr3,di3
      real :: taur,taui
!-----------------
      DATA TAUR,TAUI /-.5,.866025403784439/

      do K=1,L1
         TR2 = CC(Ido,2,K)+CC(IDO,2,K)
         CR2 = CC(1,1,K)+TAUR*TR2
         CH(1,K,1) = CC(1,1,K)+TR2
         CI3 = TAUI*(CC(1,3,K)+CC(1,3,K))
         CH(1,K,2) = CR2-CI3
         CH(1,K,3) = CR2+CI3
      end do
      if (Ido .EQ. 1) return
      IDP2 = Ido+2
      do K=1,L1
         do I=3,IDO,2
            IC = IDP2-I
            TR2 = CC(I-1,3,K)+CC(IC-1,2,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,3,K)-CC(IC,2,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,3,K)-CC(IC-1,2,K))
            CI3 = TAUI*(CC(I,3,K)+CC(IC,2,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2
            CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2
            CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3
            CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3
         end do
      end do
      return
      end

!--------------------------------------------
      subroutine R4RADB4 (Ido,L1,CC,CH,WA1,WA2,WA3)
use iso_c_binding, only: fp => c_double
      implicit none

      integer :: Ido,L1
      REAL       CC(Ido,4,L1)           ,CH(IDO,L1,4)           , &
                      WA1(*)     ,WA2(*)     ,WA3(*)

!-----------------
      integer :: i,k,ic,idp2
      real :: sqrt2,tr1,tr2,tr3,tr4,ti1,ti2,ti3,ti4
      real :: cr2,cr3,cr4,ci2,ci3,ci4
!-----------------
      DATA SQRT2 /1.414213562373095/

      do K=1,L1
         TR1 = CC(1,1,K)-CC(Ido,4,K)
         TR2 = CC(1,1,K)+CC(Ido,4,K)
         TR3 = CC(Ido,2,K)+CC(IDO,2,K)
         TR4 = CC(1,3,K)+CC(1,3,K)
         CH(1,K,1) = TR2+TR3
         CH(1,K,2) = TR1-TR4
         CH(1,K,3) = TR2-TR3
         CH(1,K,4) = TR1+TR4
      end do
      if (Ido-2.lt.0) return
      if (Ido-2.gt.0) then
         IDP2 = Ido+2
         do K=1,L1
            do I=3,IDO,2
               IC = IDP2-I
               TI1 = CC(I,1,K)+CC(IC,4,K)
               TI2 = CC(I,1,K)-CC(IC,4,K)
               TI3 = CC(I,3,K)-CC(IC,2,K)
               TR4 = CC(I,3,K)+CC(IC,2,K)
               TR1 = CC(I-1,1,K)-CC(IC-1,4,K)
               TR2 = CC(I-1,1,K)+CC(IC-1,4,K)
               TI4 = CC(I-1,3,K)-CC(IC-1,2,K)
               TR3 = CC(I-1,3,K)+CC(IC-1,2,K)
               CH(I-1,K,1) = TR2+TR3
               CR3 = TR2-TR3
               CH(I,K,1) = TI2+TI3
               CI3 = TI2-TI3
               CR2 = TR1-TR4
               CR4 = TR1+TR4
               CI2 = TI1+TI4
               CI4 = TI1-TI4
               CH(I-1,K,2) = WA1(I-2)*CR2-WA1(I-1)*CI2
               CH(I,K,2) = WA1(I-2)*CI2+WA1(I-1)*CR2
               CH(I-1,K,3) = WA2(I-2)*CR3-WA2(I-1)*CI3
               CH(I,K,3) = WA2(I-2)*CI3+WA2(I-1)*CR3
               CH(I-1,K,4) = WA3(I-2)*CR4-WA3(I-1)*CI4
               CH(I,K,4) = WA3(I-2)*CI4+WA3(I-1)*CR4
            end do
         end do
         if (MOD(Ido,2) .EQ. 1) return
      end if
      do K=1,L1
         TI1 = CC(1,2,K)+CC(1,4,K)
         TI2 = CC(1,4,K)-CC(1,2,K)
         TR1 = CC(Ido,1,K)-CC(IDO,3,K)
         TR2 = CC(Ido,1,K)+CC(IDO,3,K)
         CH(Ido,K,1) = TR2+TR2
         CH(Ido,K,2) = SQRT2*(TR1-TI1)
         CH(Ido,K,3) = TI2+TI2
         CH(Ido,K,4) = -SQRT2*(TR1+TI1)
      end do
      return
      end

!--------------------------------------------
      subroutine R4RADB5 (Ido,L1,CC,CH,WA1,WA2,WA3,WA4)
use iso_c_binding, only: fp => c_double
      implicit none

      integer :: Ido,L1
      REAL       CC(Ido,5,L1)           ,CH(IDO,L1,5)           , &
                      WA1(*)     ,WA2(*)     ,WA3(*)     ,WA4(*)

!-----------------
      integer :: i,k,ic,idp2
      real :: tr2,tr3,tr4,tr5,ti2,ti3,ti4,ti5
      real :: cr2,cr3,cr4,cr5,ci2,ci3,ci4,ci5
      real :: dr2,dr3,dr4,dr5,di2,di3,di4,di5
      real :: tr11,ti11,tr12,ti12
!-----------------
      DATA TR11,TI11,TR12,TI12 /.309016994374947,.951056516295154, &
      -.809016994374947,.587785252292473/

      do K=1,L1
         TI5 = CC(1,3,K)+CC(1,3,K)
         TI4 = CC(1,5,K)+CC(1,5,K)
         TR2 = CC(Ido,2,K)+CC(IDO,2,K)
         TR3 = CC(Ido,4,K)+CC(IDO,4,K)
         CH(1,K,1) = CC(1,1,K)+TR2+TR3
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
         CI5 = TI11*TI5+TI12*TI4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2) = CR2-CI5
         CH(1,K,3) = CR3-CI4
         CH(1,K,4) = CR3+CI4
         CH(1,K,5) = CR2+CI5
      end do
      if (Ido .EQ. 1) return
      IDP2 = Ido+2
      do K=1,L1
         do I=3,IDO,2
            IC = IDP2-I
            TI5 = CC(I,3,K)+CC(IC,2,K)
            TI2 = CC(I,3,K)-CC(IC,2,K)
            TI4 = CC(I,5,K)+CC(IC,4,K)
            TI3 = CC(I,5,K)-CC(IC,4,K)
            TR5 = CC(I-1,3,K)-CC(IC-1,2,K)
            TR2 = CC(I-1,3,K)+CC(IC-1,2,K)
            TR4 = CC(I-1,5,K)-CC(IC-1,4,K)
            TR3 = CC(I-1,5,K)+CC(IC-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2
            CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2
            CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3
            CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3
            CH(I-1,K,4) = WA3(I-2)*DR4-WA3(I-1)*DI4
            CH(I,K,4) = WA3(I-2)*DI4+WA3(I-1)*DR4
            CH(I-1,K,5) = WA4(I-2)*DR5-WA4(I-1)*DI5
            CH(I,K,5) = WA4(I-2)*DI5+WA4(I-1)*DR5
         end do
      end do
      return
      end

!--------------------------------------------
      subroutine R4RADBG (Ido,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
        use iso_c_binding, only: fp => c_double
        implicit none

      integer :: Ido,L1,IP,IDL1

      REAL       CH(Ido,L1,IP)          ,CC(IDO,IP,L1)          , &
                      C1(Ido,L1,IP)          ,C2(IDL1,IP), &
                      CH2(IDL1,IP)           ,WA(1)

!-----------------
      integer :: nbd,i,k,idp2,ipp2,ipph,j,j2,jc,ic,l,lc,ik,is,idij
      real :: ar1,ar2,ai1,ai2,ar1h,ar2h,dc2,ds2
      real :: tpi,arg,dcp,dsp
!-----------------
      DATA TPI/6.28318530717959/

      ARG = TPI/FLOAT(IP)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IDP2 = Ido+2
      NBD = (Ido-1)/2
      IPP2 = IP+2
      IPPH = (IP+1)/2
      if (Ido .LT. L1) goto 103
      do K=1,L1
         do I=1,IDO
            CH(I,K,1) = CC(I,1,K)
         end do
      end do
      goto 106
  103 do I=1,IDO
         do K=1,L1
            CH(I,K,1) = CC(I,1,K)
         end do
      end do
  106 do J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         do K=1,L1
            CH(1,K,J) = CC(Ido,J2-2,K)+CC(IDO,J2-2,K)
            CH(1,K,JC) = CC(1,J2-1,K)+CC(1,J2-1,K)
         end do
      end do
      if (Ido .EQ. 1) goto 116
      if (NBD .LT. L1) goto 112
      do J=2,IPPH
         JC = IPP2-J
         do K=1,L1
            do I=3,IDO,2
               IC = IDP2-I
               CH(I-1,K,J) = CC(I-1,2*J-1,K)+CC(IC-1,2*J-2,K)
               CH(I-1,K,JC) = CC(I-1,2*J-1,K)-CC(IC-1,2*J-2,K)
               CH(I,K,J) = CC(I,2*J-1,K)-CC(IC,2*J-2,K)
               CH(I,K,JC) = CC(I,2*J-1,K)+CC(IC,2*J-2,K)
            end do
         end do
      end do
      goto 116
  112 do J=2,IPPH
         JC = IPP2-J
         do I=3,IDO,2
            IC = IDP2-I
            do K=1,L1
               CH(I-1,K,J) = CC(I-1,2*J-1,K)+CC(IC-1,2*J-2,K)
               CH(I-1,K,JC) = CC(I-1,2*J-1,K)-CC(IC-1,2*J-2,K)
               CH(I,K,J) = CC(I,2*J-1,K)-CC(IC,2*J-2,K)
               CH(I,K,JC) = CC(I,2*J-1,K)+CC(IC,2*J-2,K)
           end do
         end do
      end do
  116 AR1 = 1.
      AI1 = 0.
      do L=2,IPPH
         LC = IPP2-L
         AR1H = DCP*AR1-DSP*AI1
         AI1 = DCP*AI1+DSP*AR1
         AR1 = AR1H
         do IK=1,IDL1
            C2(IK,L) = CH2(IK,1)+AR1*CH2(IK,2)
            C2(IK,LC) = AI1*CH2(IK,IP)
         end do
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         do J=3,IPPH
            JC = IPP2-J
            AR2H = DC2*AR2-DS2*AI2
            AI2 = DC2*AI2+DS2*AR2
            AR2 = AR2H
            do IK=1,IDL1
               C2(IK,L) = C2(IK,L)+AR2*CH2(IK,J)
               C2(IK,LC) = C2(IK,LC)+AI2*CH2(IK,JC)
            end do
         end do
      end do
      do J=2,IPPH
         do IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
         end do
      end do
      do J=2,IPPH
         JC = IPP2-J
         do K=1,L1
            CH(1,K,J) = C1(1,K,J)-C1(1,K,JC)
            CH(1,K,JC) = C1(1,K,J)+C1(1,K,JC)
         end do
      end do
      if (Ido .EQ. 1) goto 132
      if (NBD .LT. L1) goto 128
      do J=2,IPPH
         JC = IPP2-J
         do K=1,L1
            do I=3,IDO,2
               CH(I-1,K,J) = C1(I-1,K,J)-C1(I,K,JC)
               CH(I-1,K,JC) = C1(I-1,K,J)+C1(I,K,JC)
               CH(I,K,J) = C1(I,K,J)+C1(I-1,K,JC)
               CH(I,K,JC) = C1(I,K,J)-C1(I-1,K,JC)
            end do
         end do
      end do
      goto 132
  128 do J=2,IPPH
         JC = IPP2-J
         do I=3,IDO,2
            do K=1,L1
               CH(I-1,K,J) = C1(I-1,K,J)-C1(I,K,JC)
               CH(I-1,K,JC) = C1(I-1,K,J)+C1(I,K,JC)
               CH(I,K,J) = C1(I,K,J)+C1(I-1,K,JC)
               CH(I,K,JC) = C1(I,K,J)-C1(I-1,K,JC)
            end do
         end do
      end do
 132  continue
      if (Ido .EQ. 1) return
      do IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
      end do
      do J=2,IP
         do K=1,L1
            C1(1,K,J) = CH(1,K,J)
         end do
      end do
      if (NBD .GT. L1) goto 139
      IS = -Ido
      do J=2,IP
         IS = IS+Ido
         IDIJ = IS
         do I=3,IDO,2
            IDIJ = IDIJ+2
            do K=1,L1
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
            end do
         end do
      end do
      goto 143
  139 IS = -Ido
      do J=2,IP
         IS = IS+Ido
         do K=1,L1
            IDIJ = IS
            do I=3,IDO,2
               IDIJ = IDIJ+2
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
            end do
         end do
      end do
  143 return
      end


      subroutine R4_RFFT_TEST
use iso_c_binding, only: fp => c_double
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!                       VERSION 4  APRIL 1985
!
!                         A TEST DRIVER FOR
!          A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE FAST FOURIER
!           TRANSFORM OF PERIODIC AND OTHER SYMMETRIC SEQUENCES
!
!                              BY
!
!                       PAUL N SWARZTRAUBER
!
!       NATIONAL CENTER FOR ATMOSPHERIC RESEARCH  BOULDER,COLORAdo 80307
!
!        WHICH IS SPONSORED BY THE NATIONAL SCIENCE FOUNDATION
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!
!             THIS PROGRAM TESTS THE PACKAGE OF FAST FOURIER
!     TRANSFORMS FOR BOTH COMPLEX AND REAL PERIODIC SEQUENCES AND
!     CERTIAN OTHER SYMMETRIC SEQUENCES THAT ARE LISTED BELOW.
!
!     1.   RFFTI     INITIALIZE  RFFTF AND RFFTB
!     2.   RFFTF     FORWARD TRANSFORM OF A REAL PERIODIC SEQUENCE
!     3.   RFFTB     BACKWARD TRANSFORM OF A REAL COEFFICIENT ARRAY
!
!
      implicit none
!
!---------------------
!
      integer ::      ND(10)
      REAL ::         X(200)     ,Y(200)     ,W(2000)    , &
                      XH(200)

      real :: sqrt2,fn,tfn,pi,dt,sum1,sum2,arg,arg1
      real :: rftf,rftb,rftfb,cf,sum
      integer :: nns,nz,n,modn,np1,nm1,i,j,k,ns2

!---------------------
      DATA ND(1),ND(2),ND(3),ND(4),ND(5),ND(6),ND(7)/120,54,49,32,4,3,2/

      SQRT2 = SQRT(2.)
      NNS = 7
      do NZ=1,NNS
         N = ND(NZ)
         MODN = MOD(N,2)
         FN = FLOAT(N)
         TFN = FN+FN
         NP1 = N+1
         NM1 = N-1
         do J=1,NP1
            X(J) = SIN(FLOAT(J)*SQRT2)
            Y(J) = X(J)
            XH(J) = X(J)
         end do
!
!     TEST subroutineS r4RFFTI,r4RFFTF AND r4RFFTB
!
         CALL r4RFFTI (N,W)
         PI = 3.14159265358979
         DT = (PI+PI)/FN
         NS2 = (N+1)/2
         if (NS2 .LT. 2) goto 104
         do K=2,NS2
            SUM1 = 0.
            SUM2 = 0.
            ARG = FLOAT(K-1)*DT
            do I=1,N
               ARG1 = FLOAT(I-1)*ARG
               SUM1 = SUM1+X(I)*COS(ARG1)
               SUM2 = SUM2+X(I)*SIN(ARG1)
            end do
            Y(2*K-2) = SUM1
            Y(2*K-1) = -SUM2
         end do
 104     SUM1 = 0.
         SUM2 = 0.
         do I=1,NM1,2
            SUM1 = SUM1+X(I)
            SUM2 = SUM2+X(I+1)
         end do
         if (MODN .EQ. 1) SUM1 = SUM1+X(N)
         Y(1) = SUM1+SUM2
         if (MODN .EQ. 0) Y(N) = SUM1-SUM2
         CALL r4RFFTF (N,X,W)
         RFTF = 0.
         do I=1,N
            RFTF = MAX(RFTF,ABS(X(I)-Y(I)))
            X(I) = XH(I)
         end do
         RFTF = RFTF/FN
         do  I=1,N
            SUM = .5*X(1)
            ARG = FLOAT(I-1)*DT
            if (NS2 .LT. 2) goto 108
            do K=2,NS2
               ARG1 = FLOAT(K-1)*ARG
               SUM = SUM+X(2*K-2)*COS(ARG1)-X(2*K-1)*SIN(ARG1)
            end do
 108        if (MODN .EQ. 0) SUM = SUM+.5*FLOAT((-1)**(I-1))*X(N)
            Y(I) = SUM+SUM
         end do
         CALL r4RFFTB (N,X,W)
         RFTB = 0.
         do  I=1,N
            RFTB = MAX(RFTB,ABS(X(I)-Y(I)))
            X(I) = XH(I)
            Y(I) = XH(I)
         end do
         CALL r4RFFTB (N,Y,W)
         CALL r4RFFTF (N,Y,W)
         CF = 1./FN
         RFTFB = 0.
         do I=1,N
            RFTFB = MAX(RFTFB,ABS(CF*Y(I)-X(I)))
         end do

         WRITE (6,1001) N,RFTF,RFTB,RFTFB
!
 1001    format(' shortened test.f: N=',i5/ &
              '    RFTF = ',1pe11.4/ &
              '    RFTB = ',1pe11.4/ &
              '    RFTFB= ',1pe11.4)
!
      end do

      return
      end
