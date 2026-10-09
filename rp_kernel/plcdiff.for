      subroutine plcdiff(f,wk,nf)
C
C  generic finite difference operator for RPLOT calculator
C
      double precision f(nf),wk(nf)
C
      double precision fmin,fmax,dfmin
C-------------------------------------
C
C  at extrema:
C
      if(nf.eq.1) then
         f(1)=0.0
         return
      endif
C
      fmin=min(f(1),f(nf))
      fmax=max(f(1),f(nf))
      wk(1)=f(2)-f(1)
      wk(nf)=f(nf)-f(nf-1)
C
C  inside:
C
      do i=2,nf-1
         fmin=min(fmin,f(i))
         fmax=max(fmax,f(i))
         wk(i)=0.5*(f(i+1)-f(i-1))
      enddo
C
      dfmin=min(1.0d-20,1.0d-16*(fmax-fmin))
C
C  copy back; patch df=0 to protect typical d(x)/d(f) type expressions
C
      do i=1,nf
         f(i)=wk(i)
         if(abs(f(i)).lt.dfmin) then
            if(f(i).lt.0.0) then
               f(i)=-dfmin
            else
               f(i)=+dfmin
            endif
         endif
      enddo
C
      return
      end
