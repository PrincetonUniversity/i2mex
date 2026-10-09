      subroutine pbinrd(lun_warn,ichan,zfile,zdata,ndata,ngot,ier)
c
c  read a portable binary file, with given name,
c  containing zdata(1:ngot) -- single precision.
c  ngot, also read from file, must be .le. ndata = size of array passed
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
      character*(*) zfile               ! filename -- existing file
c
      integer ndata                     ! max number of data points
      real zdata(ndata)                 ! the data...
c
c  ndata must be a positive integer that will fit in 3 bytes, i.e.
c  ndata must be .le. 2**24 - 1  ~ 16,000,000
c
c  output arguments:
c
      integer ngot                      ! actual number of data pts in file
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
      call genopen(ichan,zfile,'OLD','BLOCKD',idum,ier)
      if(ier.ne.0) return
c
      nbuf=0
c
      ibytes=4                          ! byte size of local buffer
      call ubufi(ichan,buf,nbuf,ilbuf,ibytes)
c
      ifmt=ilbuf(1)                     ! =4:  code for "no compression"
      if(ifmt.ne.4) then
         write(lun_warn,*) ' ?pbinrd:  unexpected format code:  ',ifmt
         ier=1
         go to 1000
      endif
c
c  get data size and check...
c
      call ubi3dcod(ilbuf,2,ngot)       ! data size, in ilbuf(2:4).
      if(ngot.gt.ndata) then
         write(lun_warn,*) ' ?pbinrd:  passed array too small:'
         write(lun_warn,*) '  passed size = ',ndata,' need = ',ngot
         ier=1
         go to 1000
      endif
c
c  read data...
c
      ibytes=4*ngot
      call ubufi(ichan,buf,nbuf,fieee,ibytes)
c
c  reassemble byte ordering
c
      call packordr(fieee,ngot,0)
c
c  vms conversion of ieee floating point data (in unix: just copy)
c
      call vmsieee(lun_warn,zdata,ngot,fieee,0)
c
 1000 continue
c
c  close file
c
      call cftclose(ichan)
c
      return
      end
