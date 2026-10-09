!******************** START FILE FLINT.FOR ; GROUP SIGS2 ******************
!================================================================
!
!  FLINT
!   INTERPOLATE ON A 2D ARRAY
!
      REAL FUNCTION FLINT(NOUT,ARRAY,X,Y,IDENT, &
         XMIN,XMAX,XLR,NX,YMIN,YMAX,YLR,NY)
!
      DIMENSION ARRAY(NX,NY)
!
! Interpolates the table ARRAY(NX,NY) to find the value at (x,y).
! If XLR>0, then the X grid is equally spaced on a logarithmic scale:
!
!   X(i) = XMIN * (exp(XLR)) ** ( (i-1)/(NX-1) )
!
!   So X(i) ranges from XMIN to XMIN*exp(XLR), and XLR=log(xmax/xmin).
!
! For XLR<=0, the X grid is equally spaced on a linear scale from
! XMIN to XMAX.
!
! For YLR>0, Y is logarithmically spaced from YMIN to YMIN*exp(YLR)
! and for YLR<=0, Y is linearly spaced from YMIN to YMAX.
!
! NOUT is the fortran unit to write error messages to.
! IDENT is a code identifying the table ARRAY in error messages.
!
      COMMON/ZFLINT/ ZX,ZY,JDENT,IPXL,IPYL,IPXH,IPYH
!
      external zcflnt                   ! block data initialization
!
!
      ZX=X
      ZY=Y
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
      if(zy .lt. ymin) then
          zzy=0.0
      else

      if(YLR.GT.0.0) THEN
        ZLOGY=ALOG(ZY/YMIN)
        ZZY=1.0+ZLOGY/YLR * (NY-1)
      else
        ZZY=1.0+(ZY-YMIN)/(YMAX-YMIN) * (NY-1)
      end if

      end if
!
      IX=ZZX
      IY=ZZY
!  CHECK FOR POINT OUT OF BOUNDS OF ARRAY
      if(IX.GE.1) goto 10
!
      IPXL=IPXL-1
      if(IPXL.GE.0) WRITE(NOUT, 9001) IDENT,X,Y
 9001 FORMAT(/////' *************************************'/ &
        '  ATTEMPT TO INTERPOLATE OUT OF BOUNDS    IDENT=',I8/ &
        '  FLINT SIGMA*V TABLE LOOKUP, STANDARD FIXUP TAKEN'/ &
        '  ARGUMENTS:  X=',1PE10.3,'  Y=',1PE10.3/ &
        ' **************************************'//////)
      IX=1
      ZZX=1.0
!
 10   continue
      if(IY.GE.1) goto 20
      IPYL=IPYL-1
      if(IPYL.GE.0) WRITE(NOUT, 9001) IDENT,X,Y
      IY=1
      ZZY=1.0
!
 20   continue
      if(IX.LT.NX) goto 30
      IX=NX-1
      ZZX=NX
!
      if(ZX.GT.XMAX) THEN
        write(nout,9001) ident,x,y
        CALL BAD_EXIT  ! DEEMED FATAL - DMC 1988
      end if
!
 30   continue
      if(IY.LT.NY) goto 40
      IY=NY-1
      ZZY=NY
!
      if(ZY.GT.YMAX) THEN
        write(nout,9001) ident,x,y
        CALL BAD_EXIT  ! DEEMED FATAL - DMC 1988
      end if
!
!---------------------------
!  OK
!
 40   continue
      IXP1=IX+1
      IYP1=IY+1
      ZZX=ZZX-IX
      ZZY=ZZY-IY
      ZF00=(1.-ZZX)*(1.-ZZY)
      ZF01=(1.-ZZX)*ZZY
      ZF10=ZZX*(1.-ZZY)
      ZF11=ZZX*ZZY
!
!  INTERPOLATE
!
      FLINT=ZF00*ARRAY(IX,IY)+ZF01*ARRAY(IX,IYP1)+ &
              ZF10*ARRAY(IXP1,IY)+ZF11*ARRAY(IXP1,IYP1)
!
      return
      	end
!
!******************** end FILE FLINT.FOR ; GROUP SIGS2 ******************
