      subroutine trinfg2(zdisk, zdir, zrunid, zabbrev,
     >   zlabel, zunits, inx, intimes, zxname,
     >   ixref, itypref, inxref, ier)
C
C--------------------------------------------------------------
C  INITIALIZE TO LOOK AT A FUNCTION OF TIME AND X FROM ANOTHER RUN.
C
C  use trinfo to get information; verify that it is a rank 2 object;
C  verify that its named X axis corresponds to an X axis in the current
C  main RUN.
C
C  arguments are as in trinfo, except:
C    inx,intimes replaces idims(*)
C    ixref is the fcn# corresponding to the x axis data in the main run
C    itypref is the X axis type code in the main run
C    inxref is the no. of X pts for the X axis type in the main run
C
C
      use cplotr_mod
C
      character*(*) zdisk,zdir
      character*(*) zrunid
      character*(*) zabbrev
C
      character*(*) zlabel
      character*(*) zunits
C
      integer inx,intimes
      character*10 zxname
      integer ixref,itypref,inxref
C
      INTEGER       IER
C
      integer idims(8)
      character*10 zxnames(8)
C
C---------------------------------------------------------
C
      zlabel=' '
      zunits=' '
      inx=0
      intimes=0
      zxname=' '
      ixref=0
      itypref=0
      inxref=0
C
      ier=0
C
      luntrm = lunzer(0)
C
      ier = -99                           ! keep tree open
      call trinfo(zdisk, zdir, zrunid, zabbrev, itypq,
     1   zlabel, zunits, irank, idims, zxnames, ier)
C
      if(ier.eq.0) then
         if(irank.ne.2) then
            ier=999
            write(luntrm,*)
     >         ' ?inirn2:  rplot code error, expected rank 2 item.'
         else
            ixref=ifind_ordr(abr,iordrr,nfxt,zxnames(1))
            if(ixref.eq.0) then
               write(luntrm,*)
     >            ' ?inirn2:  2nd run "',zabbrev,'" x axis is "',
     >            zxnames(1),'"'
               write(luntrm,*)
     >            '  x axis unknown to current primary run.'
               ier=1000
            else
c  x axis of 2ndary run is known to primary run.
               zxname=zxnames(1)
               itypref=itypr(ixref)
               inxref=nzonex(itypref)
            endif
         endif
      endif
C
C  exit now on error
C
      if(ier.ne.0) then
         write(luntrm,9901) ier
 9901    format(' %trinfg2:  trinfo error code:',i5)
         call tconnect_close(idum)
      else
C
C  OK...
C
         inx=idims(1)
         intimes=idims(2)
C
      endif
C
      return
      end
