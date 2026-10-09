      subroutine splitn_zeff_switches(izmod,ldatzef,
     >     NLZEFM,NLZFIN,NLZFI2,NLZVBR,NLZVB2,NLZFXI,NLZSIM,NLZEFA)
C
C  set Zeff data switches according to table value (izmod=nzefmod(j))
C
C  this routine with passed arguments so that code outside TRANSP can use..
C
      implicit NONE

      integer, intent(in) :: izmod   ! integer control value
      integer, intent(in) :: ldatzef ! 1d Zeff signal data address
      !  (set .gt.0 to allow NLZFIN to be set)

      logical, intent(out) :: NLZEFM ! T for Zeff from resistivity (vsur match)
      logical, intent(out) :: NLZFIN ! T for Zeff from Zeff(t) input signal
      logical, intent(out) :: NLZFI2 ! T for Zeff from Zeff(x,t) input signal
      logical, intent(out) :: NLZVBR ! T for Zeff from VB chordal data
      logical, intent(out) :: NLZVB2 ! T for Zeff from VB array data
      logical, intent(out) :: NLZFXI ! T for Zeff from single impurity density
      logical, intent(out) :: NLZSIM ! T for Zeff from impurity densities
      logical, intent(out) :: NLZEFA ! T for Zeff profile shaping
      ! (NLZEFA used with scalar input options NLZFIN, NLZVBR)

      !---------------------------------

      NLZEFM=.FALSE.
      NLZFIN=.FALSE.
      NLZFI2=.FALSE.
      NLZVBR=.FALSE.
      NLZVB2=.FALSE.
      NLZFXI=.FALSE.
      NLZSIM=.FALSE.
      NLZEFA=.FALSE.
C
      if(izmod.eq.1) then
         NLZEFM=.TRUE.
C
      else if(izmod.eq.2) then
         NLZFIN= (ldatzef.gt.0)
C
      else if(izmod.eq.3) then
         NLZFI2=.TRUE.
C
      else if(izmod.eq.4) then
         NLZVBR=.TRUE.
C
      else if(izmod.eq.5) then
         NLZVB2=.TRUE.
C
      else if(izmod.eq.6) then
         NLZFXI=.TRUE.
C
      else if(izmod.eq.7) then
         NLZFIN= (ldatzef.gt.0)
         NLZEFA=.TRUE.
C
      else if(izmod.eq.8) then
         NLZVBR=.TRUE.
         NLZEFA=.TRUE.
C
      else if(izmod.eq.9) then
         NLZVBR=.TRUE.
         NLZFXI=.TRUE.
C
      else if(izmod.eq.10) then
         NLZSIM=.TRUE.          ! zeff from multiple impurities
C
      else if(izmod.eq.11) then
         NLZVBR=.TRUE.          ! zeff from multiple impurities normalizd by VB
         NLZSIM=.TRUE.
C
      endif
C
      return
      end
