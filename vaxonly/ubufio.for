      subroutine ubufo(lun,obuf,nbuf,out,nout)
C
      implicit NONE
C
C
C  output binary byte stream, with buffering.
C
      integer bufsiz
      parameter (bufsiz=512)  ! see also ubufi, below
C
C  all parameters are input:
C
      integer lun       ! unit of open chanel for i/o
      integer nbuf      ! portion of buffer currently in use
C
      integer*1 obuf(bufsiz) ! output buffer
      integer*1 out(*)       ! output stream segment (any size)
      integer*1 eod
C
      integer nout      ! size of segment to write (> 512 OK).
C                       ! if 0, flush buffer (EOF only)
C
      integer ii,ier,ifree,io,ip,ipl
C
C  contents of obuf and nbuf will be modified by this call.
C
      data eod/-1/
C
C------------------------------------
C
      if(nout.eq.0) then
C  buffer flush (just before closing file)
         if(nbuf.gt.-1) then
            do ii=nbuf+1,bufsiz
               obuf(ii)=eod
            enddo
         endif
         call cftwrite(lun,obuf,ier)
         nbuf=-1
         go to 1000
      endif
C
C  there is data in out array
C
C  assume there is free space in the buffer also (a full buffer would
C  already have been written).
C
      ifree=bufsiz-nbuf
C
      if(ifree.gt.nout) then
C  add to buffer (no write needed yet)
         do ii=nbuf+1,nbuf+nout
            obuf(ii)=out(ii-nbuf)
         enddo
         nbuf=nbuf+nout
         go to 1000
      else
C  fill and write first bufferfull
         do ii=nbuf+1,bufsiz
            obuf(ii)=out(ii-nbuf)
         enddo
         call cftwrite(lun,obuf,ier)
         nbuf=0
C
C  write subsequent bufferfulls
         do io=ifree+1,nout,bufsiz
            ip=io+bufsiz-1
            ipl=min(nout,ip)
            nbuf=ipl-io+1
            do ii=1,nbuf
               obuf(ii)=out(io+ii-1)
            enddo
C  write now if buffer is full (will loop through again)
            if(nbuf.eq.bufsiz) then
               call cftwrite(lun,obuf,ier)
               nbuf=0
            endif
         enddo
      endif
C
 1000 continue
      return
      end
C=====================================================================
      subroutine ubufi(lun,ibuf,nbuf,in,nin)
C
C  input binary byte stream, with buffering.
C
      implicit NONE
C
      integer bufsiz
      parameter (bufsiz=512)  ! see also ubufo, above
C
C  all parameters are input except "in(*)" which is read in or
C  filled from the buffer on this call.
C
      integer lun       ! unit of open chanel for i/o
      integer nbuf      ! ptr to first unread byte in ibuf.
C
      integer*1 ibuf(bufsiz) ! input buffer
      integer*1 in(*)        ! output stream segment (any size)
C
      integer nin       ! size of segment to read (> 512 OK).
C
C  contents of ibuf and nbuf will be changed by this call.
C
      integer iavail,ii,io,ileft,jleft,ioff,ier
C----------------------------------------
C
      if(nin.eq.0) then
         write(6,9901)
 9901    format(' %ubufi -- nin=0, programming error.')
         go to 1000
      endif
C
      if(nbuf.eq.0) then
         call cftread(lun,ibuf,ier)
         nbuf=1
      endif
C
      iavail=bufsiz-nbuf+1
      if(iavail.gt.nin) then
C
C  desired input already in buffer, just copy it; something still
C  left in buffer.
C
         do ii=1,nin
            in(ii)=ibuf(nbuf+ii-1)
         enddo
         nbuf=nbuf+nin
         go to 1000
C
      else
C
C  copy buffer and do additional reads as needed
C
         do ii=1,iavail
            in(ii)=ibuf(nbuf+ii-1)
         enddo
C
         call cftread(lun,ibuf,ier)
         nbuf=1
C
C  loop over additional buffer fulls.
C
         ileft=nin-iavail
         jleft=ileft+bufsiz
         ioff=iavail-bufsiz
         do io=1,ileft,bufsiz
            ioff=ioff+bufsiz
            jleft=jleft-bufsiz
            if(jleft.lt.bufsiz) then
               do ii=1,jleft
                  in(ii+ioff)=ibuf(ii)
               enddo
               nbuf=jleft+1
               go to 1000
            else
               do ii=1,bufsiz
                  in(ii+ioff)=ibuf(ii)
               enddo
               call cftread(lun,ibuf,ier)
               nbuf=1
            endif
         enddo
C
      endif
 1000 continue
      return
      end
