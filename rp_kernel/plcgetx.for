C*************************** Start file PlcGetX.for ; Group Rplot_sub ********
C=============================================================================
C  PlcGetX     Calculate delta x for RMAJR (Major radius X-axis)
C
C  mod dmc Aug 1999 -- use zermsg instead of uermsg for messages;
C    avoid UREAD dependence.
C
 
 
      Subroutine PlcGetX( It, dxaxis, Inx, ITypMr, Ierr)
 
      use datmgr_mod
      use cplotr_mod

      Integer Inx, Ierr, ITypMr, It
      Real dXaxis(Inx)
 
      Integer Iad
      Logical FunGot
 
C	-------------------------------------------------------------
 
 
      PLTABB = 'RMAJM'
      If (.Not. FUNGOT(IFCNMR) ) Then
          Call ZERMSG(' PlcGetx: Major Radius "RMAJM" not found.')
          Ierr = 1
          Go to 9999
      End if
 
      ITypMr = Itypr(Ifcnmr)
      InxMr  = Nzonex(ItypMr)
      If (InxMr .Ne. Inx) Then
          Call ZERMSG(' PlcGetX: Inconsistant # of pts in RMAJM.')
          Ierr = 2
          Go to 9999
      End if  ! Inx
 
C	.Read in major radius axis.
      Call DMGXOT( ItypMr, Ind1, IndDum)
 
      Call DmdLoc('RMAJM', IndMr, IsizMr, Imr)
      If (Ind1 .Eq. 0 .Or. Imr .Eq. 0)  Then
          Call Zermsg(' PlcGetx: Major radius not in memory')
          Go to 9999
      End if    ! Ind1
 
      Iad = Imr + (It-1) * Inx
      dXaxis(1) = Datbuf(Iad+1) - DatBuf(Iad)
      Do 1000 Ix=2,Inx
          dXaxis(Ix) = Datbuf(Iad+Ix-1) - Datbuf(Iad+Ix-2)
          If (Datbuf(Iad+Ix-1) .LE. Datbuf(Iad+Ix-2) ) Then
      	Call ZERMSG(
     1            ' PlcGetX: RMAJM is not monotonically increasing.')
      	Ierr = Ix
      	Go to 9999
          End If  ! Datbuf
 1000 Continue   ! Ix loop
 
9999  Continue
      Return
      End
