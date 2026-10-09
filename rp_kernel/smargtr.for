      subroutine smargtr(zstr,zval,charlist,chardflt,zchar,ier)
C
C  decode "zstr" as a real number with perhaps a special character
C  (in charlist) prepended at beginning or end
C
      character*(*) zstr                ! input string to decode
      real zval                         ! output decoded value
      character*(*) charlist            ! input list of special characters
      character*1 chardflt              ! input default special character
      character*1 zchar                 ! output special character found
      integer ier                       ! output completion code, 0 = OK
C
C  examples:
C
C      call smargtr("10.0%",zval,"RA%","A",zchar,ier)
C        returns zval=10.0, zchar="%", ier=0
C
C      call smargtr("5.2e3",zval,"RA%","A",zchar,ier)
C        returns zval=5200.0, zchar="A" (default), ier=0
C
C      call smargtr("3.14X",zval,"RA%","A",zchar,ier)
C        returns zval=0.0, zchar=" " (default), ier=1
C        (zval does not decode because X is left as its not on CHARLIST).
C
C  local:
C
      character*20 ztest
      character*1 zchar1
C
C-------------------------------------
C
      ilen=len(zstr)
      do ic=1,ilen
         if((zstr(ic:ic).ne.' ').and.
     >      (zstr(ic:ic).ne.char(9)).and.
     >      (zstr(ic:ic).ne.char(0))) then
            ic1=ic
            go to 10
         endif
      enddo
C
C  all whitespace
C
      zval=0.0
      zchar=chardflt
      ier=1
      call zermsg(' ?smargtr:  SMOOTH argument is blank or null.')
      go to 100
C
 10   continue
      do ic=ilen,ic1,-1
         if((zstr(ic:ic).ne.' ').and.
     >      (zstr(ic:ic).ne.char(9)).and.
     >      (zstr(ic:ic).ne.char(0))) then
            ic2=ic
            go to 20
         endif
      enddo
      ic2=ic1
C
C  check for special character...
C
 20   continue
      zchar=chardflt
C
      zchar1=zstr(ic1:ic1)              ! at start
      call uupper(zchar1)
      itest1=index(charlist,zchar1)
      if(itest1.gt.0) then
         zchar=zchar1
         ic1=ic1+1
         go to 50
      endif
C
      zchar1=zstr(ic2:ic2)              ! at end
      call uupper(zchar1)
      itest1=index(charlist,zchar1)
      if(itest1.gt.0) then
         zchar=zchar1
         ic2=ic2-1
         go to 50
      endif
C
 50   continue
      if(ic1.gt.ic2) then
         zval=0.0
         ier=2
         call zermsg(
     >      ' ?smargtr:  SMOOTH argument contains no value:  '//zstr)
         go to 100
      endif
C
      ztest=zstr(ic1:ic2)
      call uupper(ztest)
      ier=3
      read(ztest,'(g20.0)',err=100) zval
      ier=0
C
 100  continue
      if(ier.eq.3) then
         call zermsg(
     >      ' ?smargtr: SMOOTH argument decode failed:  '//zstr)
      endif
      return
      end
