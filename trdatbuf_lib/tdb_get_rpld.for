      REAL*8 function tdb_get_rpld(d,ztime,zRi,zYi,ier)
!
! get TF ripple magnitude in form: log(B~/B)
!
      use trdatbuf_obj
      use tdbsub_uts
      IMPLICIT NONE

      type (trdatbuf) :: d
      REAL*8,intent(in) :: ztime ! time at which to interpolate
      REAL*8,intent(in) :: zRi,zYi ! (R,Y) location, cm
      integer, intent(out) :: ier ! exit code, 0=OK
!
!============
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER isym,indr,indy,i11,i12,i21,i22,ilt,inumt,it
      integer :: lunmsg_tdb
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 zr0,zrmax,zy0,zymax,zdr,zdy,zyion,zrion,zda,zdb,zdc
      REAL*8 zdd,zr11,zy11,zp1,zp2,zdbpt,zft
!============
!
!  get ripple data at indicated location; flat extrapolation
!  actually:  log(B~/B)  (natural log)
!
!
!---------------------------------------
!
      tdb_get_rpld=-50.0_R8
      ier=0
      if(min(d%LXRPL,d%LYRPL,d%LRPL).le.0) then
         write(lunmsg_tdb(0),*) 
     >        ' ?tdb_get_rpld:  null pointers for ripple data!'
         return
      endif
!
      ZR0=d%DATBUF(d%LXRPL)
      ZRMAX=d%DATBUF(d%LXRPL+d%NXRPL-1)
      ZY0=d%DATBUF(d%LYRPL)
      ZYMAX=d%DATBUF(d%LYRPL+d%NYRPL-1)
      IF(ABS(ZY0).LE.1.0E-2_R8*ABS(ZYMAX)) THEN
         ISYM=1
      ELSE
         ISYM=0
      ENDIF
!
      ZDR=d%DATBUF(d%LXRPL+1)-ZR0
      ZDY=d%DATBUF(d%LYRPL+1)-ZY0
!
      ZYION=zYi
!
!	SET ZYION=-ZYION IF ION BELOW MIDPLANE
!
      IF((ISYM.EQ.1).AND.(ZYION.LE.0)) ZYION=-ZYION
      ZYION=MAX(ZY0,MIN(ZYMAX,ZYION))
      ZRION=MAX(ZR0,MIN(ZRMAX,zRi))
!
      INDR=1+int((ZRION-ZR0)/ZDR)
      INDR=MAX(1,MIN((d%NXRPL-1),INDR))
      do 
         if(d%datbuf(d%lxrpl+indr-1).gt.ZRION) then
            indr=indr-1
         else if(d%datbuf(d%lxrpl+indr).lt.ZRION) then
            indr=indr+1
         else
            exit
         endif
      enddo

      INDY=1+int((ZYION-ZY0)/ZDY)
      INDY=MAX(1,MIN((d%NYRPL-1),INDY))
      do 
         if(d%datbuf(d%lyrpl+indy-1).gt.ZYION) then
            indy=indy-1
         else if(d%datbuf(d%lyrpl+indy).lt.ZYION) then
            indy=indy+1
         else
            exit
         endif
      enddo

      if(d%mtrpl.eq.1) then

         ilt=d%ltime2
         inumt=d%ntime2
         call tdbsub_lookup(d%datbuf(ilt:ilt+inumt-1),inumt,ztime,
     >        it,zft)

         i11=d%lrpl + (it-1) + d%ntime2*(indR-1) + 
     >        d%ntime2*d%nxrpl*(indY-1)
         I12 = I11 + d%ntime2
         I21 = I11 + d%ntime2*d%nxrpl
         I22 = I21 + d%ntime2

         ZDA = d%DATBUF(I11)*(1-ZFT) + d%DATBUF(I11+1)*ZFT
         ZDB = d%DATBUF(I12)*(1-ZFT) + d%DATBUF(I12+1)*ZFT
         ZDC = d%DATBUF(I21)*(1-ZFT) + d%DATBUF(I21+1)*ZFT
         ZDD = d%DATBUF(I22)*(1-ZFT) + d%DATBUF(I22+1)*ZFT

      else
         ! time invariant...

         I11=(INDY-1)*d%NXRPL+INDR+d%LRPL-1
         I12=I11+1
         I21=I11+d%NXRPL
         I22=I21+1
!
!  FIND LOG OF RIPPLE AT EACH CORNER OF BOX AROUND MC ION
!
         ZDA=d%DATBUF(I11)
         ZDB=d%DATBUF(I12)
         ZDC=d%DATBUF(I21)
         ZDD=d%DATBUF(I22)
      endif
!
! ------------------------------------------------------------------------
!
      ZR11=d%DATBUF(d%LXRPL-1+INDR)
      ZDR =d%DATBUF(d%LXRPL+INDR)-ZR11
      ZY11=d%DATBUF(d%LYRPL-1+INDY)
      ZDY =d%DATBUF(d%LYRPL+INDY)-ZY11
!
      ZP1=(ZRION-ZR11)/ZDR
      ZP2=(ZYION-ZY11)/ZDY
!
      ZDBPT=ZDA+ZP1*(ZDB-ZDA)+(ZDC+(ZDD-ZDC)*ZP1)*ZP2-(ZDA+ZP1
     <   *(ZDB-ZDA))*ZP2
!
      tdb_get_rpld = ZDBPT
!
      return
      end
