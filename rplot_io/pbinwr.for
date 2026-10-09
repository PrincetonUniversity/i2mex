      subroutine pbinwr(lun_warn,ichan,zfile,zdata,ndata,ier)
c
c  write a portable binary file, with given name,
c  containing zdata(1:ndata) -- single precision.
c
      implicit NONE
c
c  input arguments:
c
      integer lun_warn                  ! fortran lun for messages
c
      integer ichan                     ! i/o channel number, btw 0 and 100
c
c  note ichan is not a fortran lun -- the underlying i/o code is in C.
c
      character*(*) zfile               ! filename -- new file
c
      integer ndata                     ! number of data points
      real zdata(ndata)                 ! the data...
c
c  ndata must be a positive integer that will fit in 3 bytes, i.e.
c  ndata must be .le. 2**24 - 1  ~ 16,000,000
c
c  output argument:
c
      integer ier                       ! completion code: 0 = normal
c--------------------------------------
c
      integer idum,nbuf,ibytes,ifmt
c
c
      integer*1 buf(512)             ! i/o buffer
      integer*1 ilbuf(4)             ! local buffer
c
      real fieee(ndata)                 ! floating data converted to ieee
c
c--------------------------------------
c
      idum=256*256*256 - 1
      if(ndata.gt.idum) then
         write(lun_warn,*) ' ?pbinwr:  data size limit exceeded.'
         write(lun_warn,*) '  ndata = ',ndata,' size limit = ',idum
         ier=1
         return
      endif
c
      call genopen(ichan,zfile,'NEW','BLOCKD',idum,ier)
      if(ier.ne.0) return
c
      nbuf=0
c
      ibytes=4                          ! byte size of local buffer
      ifmt=4                            ! =4:  code for "no compression"
      ilbuf(1)=ifmt
      call ubi3ncod(ndata,ilbuf,2)      ! data size, in ilbuf(2:4).
c
c  transfer to output buffer
c
      call ubufo(ichan,buf,nbuf,ilbuf,ibytes)
c
c  ieee conversion of floating point data (in unix: just copy)
c
      call vmsieee(lun_warn,zdata,ndata,fieee,1)
c
c  standardize byte order for write
c
      call packordr(fieee,ndata,1)
c
c  write via buffer...
c
      ibytes=4*ndata
      call ubufo(ichan,buf,nbuf,fieee,ibytes)
c
c  flush buffer
c
      call ubufo(ichan,buf,nbuf,fieee,0)
c
c  close file
c
      call cftclose(ichan)
c
      return
      end
