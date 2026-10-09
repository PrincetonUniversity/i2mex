!******************** START FILE FLIN1.FOR ; GROUP SIGS2 ******************
!================================================================
!
!  FLIN1
!   INTERPOLATE ON A 1D ARRAY
!
      REAL FUNCTION FLIN1(NOUT,ARRAY,X,IDENT,XMIN,XMAX,XLR,NX)
!
! Interpolates the table ARRAY(NX) to find the value at x.
! If XLR>0, then the X grid is equally spaced on a logarithmic scale:
!
!   X(i) = XMIN * (exp(XLR)) ** ( (i-1)/(NX-1) )
!
!   So X(i) ranges from XMIN to XMIN*exp(XLR), and XLR=log(xmax/xmin).
!
! For XLR<=0, the X grid is equally spaced on a linear scale from
! XMIN to XMAX.
!
! NOUT is the fortran unit to write error messages to.
! IDENT is a code identifying the table ARRAY in error messages.
!

!
      DIMENSION ARRAY(NX)
!
      COMMON/ZFLINT/ ZX,ZY,JDENT,IPXL,IPYL,IPXH,IPYH
!
      external zcflnt                   ! block data initialization
!
!  INTERPOLATE TO X  LOGARITHMIC SPACING IN X
!     LINEAR if XLR .LE. 0.0
!
!
      ZX=X
!
        if(zx .lt. xmin) then
            zzx=0.0
        else

      if(XLR.GT.0.0) THEN
        ZLOGX=ALOG(ZX/XMIN)
        ZZX=1.0+ZLOGX/XLR * (NX-1)
      else
        ZZX=1.0+(ZX-XMIN)/(XMAX-XMIN) * (NX-1)
      end if

      end if
!
      IX=ZZX
!
!  CHECK FOR POINT OUT OF BOUNDS OF ARRAY
      if(IX.GE.1) goto 10
!
      IPXL=IPXL-1
      if(IPXL.GE.0) WRITE(NOUT, 9001) IDENT,X
 9001 FORMAT(///' *************************************'/ &
        '  ATTEMPT TO INTERPOLATE OUT OF BOUNDS    IDENT=',I8/ &
        '  FLIN1 SIGMA*V TABLE LOOKUP, STANDARD FIXUP TAKEN'/ &
        '  ARGUMENT:  X=',1PE10.3/ &
        ' **************************************'///)
      IX=1
      ZZX=1.0
!
 10   continue
!
      if(IX.LT.NX) goto 30
      IX=NX-1
      ZZX=NX
!
      if(ZX.GT.XMAX) THEN
        write(nout,9001) ident,x
        call bad_exit
      end if
!
!---------------------------
!  OK
!
 30   continue
!
      IXP1=IX+1
      ZZX=ZZX-IX
!
!
!  INTERPOLATE
!
      FLIN1=(1.0-ZZX)*ARRAY(IX)+ZZX*ARRAY(IXP1)
!
      return
      	end
!
!******************** end FILE FLIN1.FOR ; GROUP SIGS2 ******************
!-------------------------
!  FLINT LOCAL PERM. MEMORY
!   (CF SIGS2.FOR)
!
      BLOCK DATA ZCFLNT
!
      COMMON/ZFLINT/ ZX,ZY,JDENT,IPXL,IPYL,IPXH,IPYH
!
      DATA ZX/0.0/
      DATA ZY/0.0/
      DATA JDENT/0/
!
      DATA IPXL/1/
      DATA IPYL/1/
!
      DATA IPXH/10/
      DATA IPYH/10/
!
      end
