C----------------------------------------------------------------------
C PLCCOP  REAL FUNCTION, RETURN VALUE OF OPERAND AT SPECIFIED TIME/LOCN
 
 
      REAL FUNCTION PLCCOP(Istat,IOPND,IT,IX,INX,Iwrk2, JtypMR,
     >                       iptacc)
 
 
C	Updates:
C        6/28/95 tbt Made Jtypmr an argument.
C        6/23/95 TBT Changed major radius type to variable rather
C                    than 5 (should be 7)
C        4/02/93 TbT Added Istat argument.
C                    If Istat=major radius, & IKOP = zone cen/bnd
C                    expand Ikop automatically to major radius.
C	11/23/92 TBT Added Scalar f(t) capability for $.
C
C INPUT
C  Istat - Type of non-temporal coordinate for whole calculation =
C          Maximum value of all operands.
C          -2 = constant, -1 = function of time only,
C          1 = zone center, 2= zone boundary, 7= major radius...
C  IOPND-- OPERAND NUMBER (CF RPLOT CALCULATOR COMMON)
C  IT----- TIME INDEX (PROFILE FCN TIMEBASE)
C  IX----- SPATIAL OR OTHER X AXIS INDEX
C  INX---- NUMBER OF PTS IN X AXIS INDEX
C  Iwrk2 - Address of beginning of $ in datbuf.
C  IAdl--- ADDRESS OF $(IT,IX) IN CASE NEEDED
C  JtypMr- Type for Major radius.
C  iptacc- pointer to scalar accumulator data
C
C OUTPUT
C  PLCCOP  VALUE OF OPERAND AT SPECIFIED TIME,LOCATION
C
C  OPERAND MAY BE A CONSTANT OR A FUNCTION OF TIME ONLY
C
C-----------------------------
C
      use datmgr_mod
      use cplotr_mod
      use rpcalc_mod
C
      EXTERNAL XIDENT
 
CCC	Integer IkopMR /5/     ! Type of major radius x axis.
 
C------------------------------
      Iadl = (Iwrk2-1) + (It-1)*Inx + Ix
 
C  OPERAND CLASS...
      IKOP=IKOPND(IOPND)
C
      IF(IKOP.EQ.-2) THEN
C  CONSTANT
         ZVAL=ZVOPND(IOPND)
C
      ELSE IF(IKOP.EQ.-1) THEN
C  FUNCTION OF TIME ONLY; INTERPOLATION may be NEEDED TO PROFILE TIMEBASE
         CALL DMDLOC('%F(T)',JTIM,ISIZE,IPT)
         IFCN=IAOPND(IOPND)
         if(istat.lt.0) then
c  output is also scalar, also on scalar timebase (no interpolation)
            IF(IFCN.EQ.0) THEN
               ZVAL=DATBUF(iptacc+it-1)
            ELSE
               ZVAL=DATBUF(IPT+(IFCN-1)*NTT+IT-1)
            ENDIF
         else
c  output is on profile timebase, interpolation needed.
            ZT=TIME3(IT)
            CALL XINTER(XIDENT,ZT,TIME,NTT,
     >         IT0,IT0P1,ZTI,ZTIC,IEX)
            if(ifcn.eq.0) then
               iat0=iptacc+it0-1
               iat0p1=iptacc+it0p1-1
            else
               IAT0=IPT+(IFCN-1)*NTT+IT0-1
               IAT0P1=IPT+(IFCN-1)*NTT+IT0P1-1
            endif
            ZVAL=DATBUF(IAT0)*ZTIC+DATBUF(IAT0P1)*ZTI
         End If                         ! Ifcn
C
      ELSE If (Ikop .Eq. Istat  .Or.
     1      (Ikop .Eq. 1 .And. Istat .Eq. 2) .Or.
     2      (Ikop .Eq. 2 .And. Istat .Eq. 1)) Then
C     .Compatible PROFILE FUNCTION
         IFCN=IAOPND(IOPND)
         IF(IFCN.EQ.0) THEN
            ZVAL=DATBUF(IAdl)           ! Handle accumulator, $
         ELSE
            ILOC=LOCD(IFCN)+(IT-1)*INX+IX-1
            ZVAL=DATBUF(ILOC)
         ENDIF
C
CCC	Else If (Istat.Eq.IkopMr .And. Ikop.Eq.2) Then
      Else If (Istat.Eq.JtypMr .And. Ikop.Eq.2) Then
C	    .Handle case of expanding zone bnd profile to major radius.
 
C	  Ikopnd(Iopnd) = Ikopmr     ! Changing to major radius - do in PLCFXT.
         IFCN=IAOPND(IOPND)
         Inxzone  = NzoneX(Ikop)
         Ixzone   = Abs(Ix-(Inx/2)-1)
         IF(IFCN.EQ.0) THEN
C           .handle accumulator
            If (Ixzone .Gt. 0) Then
C              .Map Zone bnd accumulator value to major radius.
               ILOC=Iwrk2 +(IT-1)*INXzone+IXzone-1
               ZVAL=DATBUF(ILOC)
            Else
C	       .Calculate center point in major radius axis.
               ILOC=Iwrk2 +(IT-1)*INXzone
               ZVal=Datbuf(Iloc)+ .333333*(Datbuf(Iloc)-Datbuf(Iloc+1))
            End If                      ! Ixzone
 
         ELSE                           ! Ifcn .Ne. 0
C	    . Not accumulator.
            If (Ixzone .Gt. 0) Then
C              .Map Zone bnd value to major radius.
               ILOC=LOCD(IFCN)+(IT-1)*INXzone+IXzone-1
               ZVAL=DATBUF(ILOC)
            Else
C	       .Calculate center point in major radius axis.
               ILOC=LOCD(IFCN)+(IT-1)*INXzone
               ZVal=Datbuf(Iloc)+ .333333*(Datbuf(Iloc)-Datbuf(Iloc+1))
            End If                      ! Ixzone
         ENDIF                          ! Ifcn
 
CCC	Else If (Istat.Eq.5 .And. Ikop.Eq.1) Then   ! tbt 6/95
      Else If (Istat.Eq.JtypMr .And. Ikop.Eq.1) Then
C	    .Handle case of expanding zone centered profile to major radius.
C           NOTE: Values are interpolated to zone boundaries for major radius
C                 conversion. See subroutine XIntzB.
 
C	  Ikopnd(Iopnd) = Ikopmr      ! Changing to major radius - Do in PLCFXT
         IFCN=IAOPND(IOPND)
         Inxzone  = NzoneX(Ikop)
         Ixzone   = Abs(Ix-(Inx/2)-1)
         IF(IFCN.EQ.0) THEN
C           .handle accumulator
            If (Ixzone .Gt. 0) Then
C              .Map Zone centered accumulator value to major radius.
               ILOC=Iwrk2 +(IT-1)*INXzone+IXzone-1
C	       ZVAL=DATBUF(ILOC)
C              .Convert from zone center to zone boundaries.
               If (Ixzone .Lt. Inxzone) Then
                  Zval = 0.5 * (Datbuf(Iloc)+Datbuf(Iloc+1))
               Else
C     Extrapolate at edge.
                  Zval = Datbuf(Iloc)+0.5*(Datbuf(Iloc)-Datbuf(Iloc-1))
               End If                   ! Ixzone
            Else
C	       .Calculate center point in major radius axis.
               ILOC=Iwrk2 +(IT-1)*INXzone
C              .Convert from zone center to zone boundaries.
C	       ZVal=Datbuf(Iloc)+ .333333*(Datbuf(Iloc)-Datbuf(Iloc+1))
               Zval = 0.5 * (Datbuf(Iloc)+Datbuf(Iloc+1)) +
     1            0.333333 * 0.5 * (Datbuf(Iloc)-Datbuf(Iloc+2))
            End If                      ! Ixzone
 
         ELSE                           ! Ifcn .Ne. 0
C	    . Not accumulator.
            If (Ixzone .Gt. 0) Then
C              .Map Zone bnd value to major radius.
               ILOC=LOCD(IFCN)+(IT-1)*INXzone+IXzone-1
C	       ZVAL=DATBUF(ILOC)
C              .Convert from zone center to zone boundaries.
               If (Ixzone .Lt. Inxzone) Then
                  Zval = 0.5 * (Datbuf(Iloc)+Datbuf(Iloc+1))
               Else
C		  Extrapolate at edge.
                  Zval = Datbuf(Iloc)+0.5*(Datbuf(Iloc)-Datbuf(Iloc-1))
               End If                   ! Ixzone
            Else
C	       .Calculate center point in major radius axis.
               ILOC=LOCD(IFCN)+(IT-1)*INXzone
C	       ZVal=Datbuf(Iloc)+ .333333*(Datbuf(Iloc)-Datbuf(Iloc+1))
               Zval = 0.5 * (Datbuf(Iloc)+Datbuf(Iloc+1)) +
     1            0.333333 * 0.5 * (Datbuf(Iloc)-Datbuf(Iloc+2))
            End If                      ! Ixzone
         ENDIF                          ! Ifcn
 
      Else
C	    .It should be impossible to get here, but just in case.....
         Call zermsg('??Plccop: Illegal function type', Istat, 0)
         Call abortt
         Zval = 0.
      ENDIF
C
      PLCCOP=ZVAL
      RETURN
      END
