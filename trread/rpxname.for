      subroutine rpxname(i,zxname)
C
C  fetch x axis function name, given x axis type code
C
      use cplotr_mod

      integer i                         ! x axis type code (in)
      character*10 zxname               ! x axis function name (out)
C
C  if an invalid type code is given, zxname=' ' is returned.
C
C-----------------------------
C
      zxname=' '
      if((1.le.i).and.(i.le.nxr)) then
         zxname=xndabb(i)
      endif
C
      return
      end
