      subroutine plcgarg(inumarg,zstr)
      use rpcalc_mod
C
C  RPLOT calculator command utility routine:  get command argument
C
      integer inumarg                   ! input:  argument number
      character*(*) zstr                ! output: value (could be default)
C
C------------------------------
C
      zstr='error'
C
      if(kcmd.gt.0) then
         if((inumarg.ge.1).and.(inumarg.le.ncmdargs(kcmd))) then
            if(kposarg(inumarg).gt.0) then
               ik1=kposarg(inumarg)
               ik2=ik1+klenarg(inumarg)-1
               zstr=cmdbuf(ik1:ik2)
            else
               zstr=rppadfs(inumarg,kcmd)
            endif
         endif
      endif
C
      return
      end
 
