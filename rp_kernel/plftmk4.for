C----------------------------------------------------------------------
C  INSERT SCALAR FCN IN RPLOT DATA BUFFER FOR LATER REFERENCE
C
 
      SUBROUTINE PLFTMK4(ZTIME,ZDATA,INPTS,ZLBL,ZUNS, ZAbr)
 
      use datmgr_mod
      use cplotr_mod

C  update dmc 3 Aug 1999 -- code common to PLFTMK*.FOR placed into
C  subroutines.
C
C	Created from PLFTMK, 12/3/92 by TBT. This version does not
C       request information from the screen. All information to form
C       the new scalar is in the arguments. Created to load TIME
C       into the scalar data base of RPLOT.
C	Called from RPLOT.for
C
C  this routine assumes there is space for the data in the f(t) buffer;
C  no check is performed.  Also the abbreviation is assumed to be valid.
C
      Character*(*) ZAbr   ! Input scalar name abbreviation.
      CHARACTER*(*) ZLBL   ! Input label
      CHARACTER*(*) ZUNS   ! Input units
      REAL ZTIME(INPTS)  ! INPUT TIME SEQUENCE
      REAL ZDATA(INPTS)  ! INPUT TIME SERIES DATA
C
      CHARACTER*21 ZINPUT
C
      REAL ZANSI(60)
C
      EXTERNAL XIDENT
C
C--------------------------------
C
      LUNT=lunzer(0)
C
C  CHECK THAT THE INPUT TIME DATA IS IN ORDER
C
      call plftmck('PLFTMK4',LUNT,ztime,zdata,inpts,IORDR,ZAMAX)
C
      ZInput = ZAbr
C
      call plftmake(lunt,zinput,zlbl,zuns,iordr,zamax,
     >   ztime,zdata,inpts)
C
      return
      end
