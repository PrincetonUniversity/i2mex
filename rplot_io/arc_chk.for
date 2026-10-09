      subroutine arc_chk(runid,iwait)
c
c  *avoiding* ureadsub, ask if user will wait for a run to be restored
c
c    return iwait=1 to wait
c                =0 not to wait
c
      character*(*) runid
      integer iwait
c
c-----------------------------------------
c  local:
c
      integer ibatch                    ! =1:  batch mode
c
      character*1 zans
c-----------------------------------------
c
      call isbatch(ibatch)
      if(ibatch.eq.1) then
c
c  no terminal attached...
c
         write(6,*)
     >      ' %arc_chk:  '//runid//' is off-line, waiting...'
         iwait=1
c
      else
c
         write(6,*)
     >      ' %arc_chk:  '//runid//' is off-line.'
         write(6,*) ' ---> wait for run to be restored? (Y/N):'
         read(5,'(A1)') zans
         if((zans.eq.'y').or.(zans.eq.'Y')) then
            write(6,*) ' ...waiting... '
            iwait=1
         else
            write(6,*) ' ...restore was requested, try again later.'
            iwait=0
         endif
c
      endif
c
      return
      end
