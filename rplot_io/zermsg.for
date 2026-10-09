      subroutine zermsg(zmsg)
      use rpcalc_mod
C
      character*(*) zmsg
C
      write(lunrpc,'(1x,A)') zmsg
C
      return
      end
C--------------------------------
      integer function lunzer(idum)
C
      use rpcalc_mod
C
      lunzer=lunrpc
C
      return
      end
