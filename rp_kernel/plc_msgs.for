      subroutine plc_msgs(ilun,zfile)

      use datmgr_mod
      use rpcalc_mod
C
C  redirect rplot calculator messages to file 'zfile' on i/o unit 'ilun'.
C
C  to have messages come out on stdout (which is the default, cf
C  rpcaldat.for), call plc_msgs(6,' ').
C
      integer ilun                      ! i/o unit number to use
      character*(*) zfile               ! file to write
C
C----------------------------------------------
C  verify initialization...
C
      call initcpl
C
      lunrpc=ilun
      lundmo=ilun
      if(zfile.ne.' ') then
         call genopen(lunrpc,zfile,'UNKNOWN','ASCII',idum,ier)
         if(ier.ne.0) then
            lunrpc=6
            lundmo=6
            write(6,1001) zfile
 1001       format(' ?plc_msgs -- could not open for write:  ',a/
     >         '  rplot calculator messages written to stdout lun=6')
         endif
      endif
C
      return
      end
 
 
