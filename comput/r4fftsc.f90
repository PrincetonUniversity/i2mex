subroutine r4fftsc  (a,n,st,ct,ifail)
  use iso_c_binding, only: fp => c_double
!-----------------------------------------------------------------------
!  DMC -- real precision IMSL replacement routine
!         locally allocated, saved workspace objects for each problem
!         size (n) that is seen.  This may not be thread safe;
!         this is only intended to replace copied IMSL routines in some
!         serial TRANSP-related data processing applications.  For codes
!         requiring high performance FFT, use FFTW 2.1.5.
!
!--------------------------------
!   some IMSL comments retained...
!--------------------------------
!
!   purpose             - compute the sine and cosine transforms of
!                           a real valued sequence
!
!   usage               - call fftsc (a,n,st,ct,iwk,wk,cwk)
!
!   arguments    a      - input real vector of length n which
!                           contains the data to be transformed.
!                n      - input number of data points to be transformed.
!                           n must be a positive even integer ::.
!                st     - output real vector of length n/2+1
!                           containing the coefficients of the
!                           sine transform.
!                ct     - output real vector of length n/2+1
!                           containing the coefficients of the
!                           cosine transform.
!
!              ifail    - 0 on exit: normal
!                       - 1 on exit: could not allocate workspace
!
!   notation            - information on special notation and
!                           conventions is available in the manual
!                           introduction or through imsl routine uhelp
!
!   remarks  1.  fftsc computes the sine transform, st, according
!                to the following formula;
!
!                  st(k+1) = 2.0 * sum from j = 0 to n-1 of
!                            a(j+1)*sin(2.0*pi*j*k/n)
!                  for k=0,1,...,n/2 and pi=3.1415...
!
!                fftsc computes the cosine transform, ct, according
!                to the following formula;
!
!                  ct(k+1) = 2.0 * sum from j = 0 to n-1 of
!                            a(j+1)*cos(2.0*pi*j*k/n)
!                  for k=0,1,...,n/2 and pi=3.1415...
!            2.  the following relationship exists between the data
!                and the coefficients of the sine and cosine transform
!
!                  a(j+1) = ct(1)/(2*n) + ct(n/2+1)/(2*n)*(-1)**j +
!                           sum from k = 1 to n/2-1 of
!                            (ct(k+1)/n*cos((2.0*pi*j*k)/n) +
!                             st(k+1)/n*sin((2.0*pi*j*k)/n))
!                  for j=0,1,...,n-1 and pi=3.1415...
!
!
!-----------------------------------------------------------------------
!
      implicit none
!
!                                  specifications for arguments
      integer ::            n
      real               a(n),st(*),ct(*)
      integer ::            ifail
!
!-----------------------------------------------------------------------
!
      type :: fftwk
         integer :: n_size
         real, dimension(:), pointer :: wka => NULL()
      end type fftwk
!
      real, dimension(:), pointer :: wka => NULL()
!
      integer, SAVE :: n_save = 0
!
!  collection of work arrays:
!
      integer, parameter :: nmax=50
      type (fftwk), dimension(nmax), SAVE :: fw
!
!-----------------------------------------------------------------------
      integer :: isize,in2p1,ino2,k1,k2,i,imatch
      real, dimension(:), allocatable :: awk
!-----------------------------------------------------------------------
!
      ifail = 0
!
      if(n.le.0) return
!
      in2p1 = n/2 + 1

      st(1:in2p1) = 0.0
      ct(1:in2p1) = 0.0
!
!-------------------------------
!
      imatch=0
      do i=1,n_save
         if(fw(i)%n_size.eq.n) then
            imatch=i
            exit
         end if
      end do

      if(imatch.eq.0) then
         if(n_save.eq.nmax) then
            write(6,*) ' ?r4fftsc: nmax exceeded. '
            ifail=1
            return
         end if

         n_save = n_save + 1
         fw(n_save)%n_size = n

         isize = 2*n + 30
         allocate(fw(n_save)%wka(isize)); fw(n_save)%wka = 0.0

         call r4rffti(n,fw(n_save)%wka)

         imatch = n_save
      end if

      wka => fw(imatch)%wka
!
!-------------------------------
!
      allocate(awk(n))
      awk = a
!
      call r4rfftf(n,awk,wka)
!
      ct(1)=awk(1)*2.0
      st(1)=0.0
!
      k1=0
      k2=1
!
      ino2 = (n+1)/2
!
      do i=2,ino2
         k1=k1+2
         k2=k2+2

         ct(i) = awk(k1)*2.0
         st(i) = -awk(k2)*2.0
      end do
!
      if(ino2.lt.in2p1) then
         k1=k1+2
         ct(i)=awk(k1)*2.0
      end if

      return
      end
