C-----------------------------------------------------------------------
C  TDB_UNMAP -- INTERPOLATE THE UNMAPPED INPUT PROFILE DATA
C
      SUBROUTINE TDB_UNMAP(d,t, INRIA,IHECE,ZTIME,ILX,INX,ILF, ierr)
C
C  THIS SUBROUTINE IS CALLED FROM PROFLI -- THE PASSED ARGUMENTS
C  ARE ALL PROFLI PASSED ARGUMENTS, DESCRIBED IN THE COMMENTS NEAR
C  THE TOP OF PROFLI.FOR
C
C  INPUT:  INRIA=IABS(INRI),IHECE,ZTIME,ILX,INX,ILF
C  OUTPUT:  t%datusym(...)
C
      use trdatbuf_obj
      use trdatbuf_aux
      IMPLICIT NONE

      type (trdatbuf) :: d
      type (profget) :: t

      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
C
!============
! idecl:  explicitize implicit INTEGER declarations:
      INTEGER ihece,ilx,inx,ilf,inria,ierr,iixtrp,icenr,i,ii,irefl
      INTEGER icen,inc,incf,ippin,ippout,j,ipp,i3,i2,i1,imin,igap
      INTEGER insave,ifreq,ipatch,iisave
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 ztime,zxi,zsign,zpf,zslp,zfbin,zfbout,zfcur,zfnew,zdr
      REAL*8 zr1,zd1,zr2,zd2,zrden
!============
C
C  INRIA SPECIFIES THE DATA MAPPING OPTION
C
C  IHECE.EQ.0 IF THIS IS NORMAL PROFILE DATA
C  IHECE.EQ.<harmonic> IF THIS IS ECE DATA
C
C  ZTIME -- INTERPOLATION TARGET TIME
C
C  ILX,INX,ILF -- POINTERS TO THE DATA
C
C  LOCAL STUFF:
      real*8, dimension(:), allocatable :: zxl,ztemp
      real*8 :: rminrd,rmajrd
C
      LOGICAL ILDBUG  ! .TRUE. FOR EXTRA DEBUG OUTPUT
      DATA ILDBUG/.FALSE./
C
      integer isave(2,10),nrsym,nrmaj,nzones,izp1
      integer :: nonlin, lunmsg_tdb
C
C-----------------------------------------------------------------------
C
      IERR=0
      IIXTRP=0
C
      t%ecegap = 0
C
      nonlin = lunmsg_tdb(0)
C
      nrsym = t%nrsym
      nrmaj = t%nrmaj
      nzones = t%nzones
      izp1=nzones+1
      rminrd=0.5_R8*(t%rmajmp(nrmaj)-t%rmajmp(1))
      rmajrd=0.5_R8*(t%rmajmp(nrmaj)+t%rmajmp(1))
C
      allocate(zxl(nrsym),ztemp(nrsym))
C
      IF(IHECE.EQ.0) THEN
C  NON-ECE DATA
         ICENR=2*NZONES+3
         DO I=1,NRSYM
            IF((INRIA.EQ.5).OR.(INRIA.EQ.8)) THEN
               ZXL(I)=t%XIRSYM(I)
               if(inria.eq.8) zxl(i)=abs(zxl(i))*zxl(i)
            else if((inria.eq.6).or.(inria.eq.7)) then
               zxi=t%xirsym(i)
               if(zxi.lt.0.0_R8) then
                  zsign=-1.0_R8
                  zxi=-zxi
               else
                  zsign=+1.0_R8
               endif
               if(zxi.ge.1.0_R8) then
                  zpf=t%plflxg(izp1)+(t%plflxg(izp1)-t%plflxg(nzones))*
     1               (zxi-1.0_R8)/(t%xibdys(izp1)-t%xibdys(nzones))
               else
                  ii=1+nzones*zxi
                  ii=min(nzones,max(1,ii))
 88               continue
                  if(zxi.lt.t%xibdys(ii)) then
                     ii=ii-1
                     go to 88
                  else if(zxi.gt.t%xibdys(ii+1)) then
                     ii=ii+1
                     go to 88
                  endif
                  zpf=t%plflxg(ii)+(t%plflxg(ii+1)-t%plflxg(ii))*
     1               (zxi-t%xibdys(ii))/(t%xibdys(ii+1)-t%xibdys(ii))
                  zpf=max(0.0_R8,zpf)
               endif
               zxl(i)=zsign*zpf/t%plflxg(izp1)
               if(inria.eq.7) zxl(i)=zsign*sqrt(zpf/t%plflxg(izp1))
            ELSE IF(INRIA.EQ.4) THEN
               IREFL=ICENR-(I-ICENR)
               ZXL(I)=0.5_R8*(t%RMJSYM(I)-t%RMJSYM(IREFL))/RMINRD
            ELSE
               ZXL(I)=(t%RMJSYM(I)-RMAJRD)/RMINRD
            ENDIF
         ENDDO
C
         ZSLP=0.0_R8
         CALL INT2D(d%DATBUF(d%LTIME2),d%NTIME2,
     1			d%DATBUF(ILX),INX,
     2			d%DATBUF(ILF),d%NTIME2,INX,
     3			ZTIME,
     4			ZXL, NRSYM,
     4			t%datusym,
     5			ZSLP,
     6			IIXTRP, ILDBUG, IERR, NONLIN)
C
      ELSE
C
C  ECE DATA WITH RADIUS TO FREQUENCY MAP
C
        ICEN=NZONES+1
        ICENR=2*NZONES+3
        ZXL(ICENR)=28.0_R8*IHECE*t%BMIDP(icen)
C
        INC=0
        INCF=0
        IPPIN=ICEN
        IPPOUT=ICEN
        DO J=2,izp1
C  BUILD DBL RESOLUTION FREQUENCY VECTOR -- INTERPOLATION TARGET
C  ORDER IS REVERSED FROM ORDER OF MAJOR RADII
          INC=INC+1
          INCF=INCF+1
C
          IPP=IPPIN
          IPPIN=ICEN-INCF
          ZFBIN=0.5_R8*(t%bmidp(IPPIN)+t%bmidp(IPP))
          ZXL(ICENR+INC)=28.0_R8*IHECE*ZFBIN
C
          IPP=IPPOUT
          IPPOUT=ICEN+INCF
          ZFBOUT=0.5_R8*(t%bmidp(IPPOUT)+t%bmidp(IPP))
          ZXL(ICENR-INC)=28.0_R8*IHECE*ZFBOUT
C
          INC=INC+1
C
          ZXL(ICENR+INC)=28.0_R8*IHECE*t%bmidp(ippin)
          ZXL(ICENR-INC)=28.0_R8*IHECE*t%bmidp(ippout)
C
        ENDDO
C
        ZXL(2)=ZXL(3)*T%RMJSYM(2)/T%RMJSYM(3)
        ZXL(1)=ZXL(3)*T%RMJSYM(1)/T%RMJSYM(3)
C
        I3=4*NZONES+3
        I2=4*NZONES+4
        I1=4*NZONES+5
        ZXL(I2)=ZXL(I3)*T%RMJSYM(I2)/T%RMJSYM(I3)
        ZXL(I1)=ZXL(I3)*T%RMJSYM(I1)/T%RMJSYM(I3)
C
C  DMC Oct 17 1995
C  enforce monotonicity of ZXL (ECE frequencies); print out warning
C  if non-monotonicity occurs.  This indicates a B(R) singularity.
C
        imin=1
        zfcur=zxl(1)
        igap=0
        insave=0
C
        do ifreq=2,nrsym
           zfnew=1.000001_R8*zfcur
           if(zxl(ifreq).le.zfnew) then
              igap=ifreq
              zxl(ifreq)=zfnew
              ipatch=0
           else
              ipatch=1
           endif
           if(ifreq.eq.nrsym) ipatch=1
           if(ipatch.eq.1) then
              if(igap.gt.0) then
                 zdR=t%rmjsym(nrsym+1-imin)-t%rmjsym(nrsym+1-igap)
                 write(nonlin,9900) zdR
                 t%ecegap=max(t%ecegap,zdR)
                 insave=insave+1
                 isave(1,insave)=imin+1
                 isave(2,insave)=igap
                 igap=0
              endif
              imin=ifreq
           endif
           zfcur=zxl(ifreq)
        enddo
 9900   format(' %TDB_UNMAP -- B(R) singularity, ECE gap of approx. ',
     1       1pe11.4,' cm.')
C
        if(igap.gt.0) then
           write(nonlin,9991)
           ierr=99
           go to 900
 9991      format(
     >   ' ?TDB_UNMAP -- highest frequency data pt in B(R) singularity')
        endif
C
        ZSLP=0.0_R8
        CALL INT2D(d%DATBUF(d%LTIME2),d%NTIME2,
     1			d%DATBUF(ILX),INX,
     2			d%DATBUF(ILF),d%NTIME2,INX,
     3			ZTIME,
     4			ZXL, NRSYM,
     4			ZTEMP,
     5			ZSLP,
     6			IIXTRP, ILDBUG, IERR, NONLIN)
C
C  ORDER REVERSAL AGAIN FOR OUTPUT PLOT
C
        DO I=1,NRSYM
          t%datusym(I)=ZTEMP(NRSYM+1-I)
        ENDDO
C
C  patch singularities
        do iisave=1,insave
           i1=nrsym+1-isave(2,iisave)
           i2=nrsym+1-isave(1,iisave)
           if((i1.eq.1).or.(i2.eq.nrsym)) then
              write(nonlin,9992)
 9992         format(
     >             ' ?TDB_UNMAP -- B(R) singularity at bdy of ECE data'/
     >             '           mapping depends on extrapolation.')
              if(i1.eq.1) i1=2
              if(i2.eq.nrsym) i2=nrsym-1
           endif
C
           zr1=t%rmjsym(i1-1)
           zd1=t%datusym(i1-1)
           zr2=t%rmjsym(i2+1)
           zd2=t%datusym(i2+1)
C
           zrden=1.0_R8/(zr2-zr1)
           do ii=i1,i2
              t%datusym(ii)=zd1+(t%rmjsym(ii)-zr1)*(zd2-zd1)*zrden
           enddo
        enddo
C
      ENDIF  ! ECE / not ECE
C
      IF(IERR.GT.0) THEN
        if(IHECE.gt.0) then
           write(nonlin,9901) t%item
        else
           write(nonlin,9902) t%item
        endif
 9901   format(' ?TDB_UNMAP -- INT2D error processing ECE data ',a)
 9902   format(' ?TDB_UNMAP -- INT2D error processing non-ECE data ',a)
      ENDIF
C
 900  continue
      deallocate(zxl,ztemp)
      RETURN
      END
