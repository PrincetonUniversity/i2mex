      integer function idcod_splitn(zint,ierr)
C
C  decode an integer -- splitn subroutine
C
C  input:
C
      character*(*) zint                !character string containing integer
C
C  output:
C
      integer ierr                      !exit code, 0=OK
C
C     function value:  the value of the integer.
C
      character*10 zintbuf
C-----------------------------
C
      ierr=1
      idcod_splitn=0
C
      ilen=len(zint)
C
C  first nonblank:
C
      do ic=1,ilen
         if((zint(ic:ic).ne.' ').and.(zint(ic:ic).ne.char(9))) go to 10
      enddo
C
C  all blank (this is an error)
C
      return
C
 10   continue
      ic1=ic
C
C  last nonblank
C
      do ic=ilen,1,-1
         if((zint(ic:ic).ne.' ').and.(zint(ic:ic).ne.char(9))) go to 20
      enddo
C
 20   continue
      ic2=ic
C
C  OK try to read as a number
C
      ilc=ic2-ic1+1
      if(ilc.gt.len(zintbuf)) return    ! error:  field too long
C
      ia=len(zintbuf)-ilc+1
      zintbuf=' '
      zintbuf(ia:ia+ilc-1)=zint(ic1:ic2)
C
      read(zintbuf,'(I10)',err=90) ians
      go to 100
C
C  error
C
 90   continue
      return
C
C  OK
C
 100  continue
      ierr=0
      idcod_splitn=ians
      return
      end
