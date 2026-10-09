c
c========================================================================
c
      subroutine cmprs2(bitems, nitems, pack, maxp, bit, ierr)
c
c
c dmc 27 mar 1996
c  this is a variant on cmprsi, for a portable compressed binary
c  ufile application.  The input stream is in bytes rather than I*2,
c  to avoid byte ordering dependencies.  The input stream is always
c  2 bytes (16 bits) wide -- not controllable as was the case for
c  cmprsi.  Also, "bit" is taken to be zero on entry.
c
c  the output stream is also a byte sequence -- 4 bytes wide.
c  standardize the byte ordering of the encoding here, for portability.
c
c marilee's original comments:
c
c+an algorithm to compress waveforms by expressing them
c as differences and storing only the number of bits needed
c to express the values over a power of two range.
c look ahead is limited to two levels but
c rescan is restarted with each output
c
c bitems  byte stream input, 2 x nitems bytes.
c nitems  number of value-pairs to be packed
c pack	  output packed array
c maxp	  dimension of pack
c bit	  starting left-to-right bit number in pack
c	  returns last occupied (starts at 0)
c	  assumed 0 on entry, pack is zeroed; note values are ored into pack
c ierr	  returned error code (0 = ok)
c	  15 if  unpacked value too large to store
c	  14 if  dimension of pack too small
c
c*25-mar-83 ---cmprsi.for, xpandi.for---
c*16-nov-83  m. thompson (ppl) - return error if nbits > 16 rather than stop
c	check dimension of pack
c*24-jan-84  m. thompson (ppl) - store count word in 16 bits vs 12
c	check for data values too large to store
c*18-jan-93  m. thompson (ppl) - rearrange logic in cmprsi
c*26-feb-92  m. thompson (ppl) - replace macro packi/unpki with in-line
c	code using shifts and masks; increment bit size by 2 (vs 1)
c
c*27-mar-96  dmc -- cmprs2/xpand2 variants
c
c-----------------------------------------------------------------------
c
      implicit integer(a-z)
c
      parameter (shift=16)	!was 12 (2/93, met)
      parameter (two=2*shift,ko256=shift+256)
c
      integer*1 bitems(2,nitems)
c
      integer pack(maxp)
c
      integer*4 iold,it
c
      integer power(17)
c
      integer sig(65)
c
      integer*2 is(ko256)		! stored differences
c
      integer*4 masks(0:16)
c
      data power/
     >    0,1,2,4,8,16,32,64,128,256,512,
     >    1024,2048,4096,8192,16384,32768/
      data sig/
     >    0,2,2,2*4,4*4,8*6,16*6,32*7/
#ifdef __F90
      data masks/ Z"00000000",
     >	  Z"00000001",Z"00000003",Z"00000007",Z"0000000f",
     >	  Z"0000001f",Z"0000003f",Z"0000007f",Z"000000ff",
     >	  Z"000001ff",Z"000003ff",Z"000007ff",Z"00000fff",
     >	  Z"00001fff",Z"00003fff",Z"00007fff",Z"0000ffff"/
#else
      data masks/ '00000000'x,
     >	  '00000001'x,'00000003'x,'00000007'x,'0000000f'x,
     >	  '0000001f'x,'0000003f'x,'0000007f'x,'000000ff'x,
     >	  '000001ff'x,'000003ff'x,'000007ff'x,'00000fff'x,
     >	  '00001fff'x,'00003fff'x,'00007fff'x,'0000ffff'x/
#endif
c-----------------------------------------------------------------------
c
      nb=16  ! always 16 in this variant
      bit=0  ! assumed 0 initially
c
      if(bit.eq.0) then
        do ii=1,maxp
          pack(ii)=0
        enddo
      endif
c
c first difference to remove most of value
      offset=power(nb+1)
      mask=2*offset-1
      imask=mask
      x0=0
      iold=0
      kptr=bit/32+1
      kpos=mod(bit,32)
      kpos=32-kpos
c loop 1 to nitems+1, limited to 256 items
      do while (x0.lt.nitems)
        nn = min(nitems-x0,ko256)
        y1=-1		!number of bits
        y2=-1
        y12=0
        x2=0		!count at this bit size
        do x3=0,nn
          if(x3.eq.nn) then
c hit the count wall
            y3=64
            go to 7
          endif
c significance
c  dmc -- form integer from byte pair
          bytehi=bitems(1,x3+x0+1)
          bytelo=bitems(2,x3+x0+1)
c
          if(bytehi.lt.0) bytehi=bytehi+256
          if(bytelo.lt.0) bytelo=bytelo+256
          it=256*bytehi+bytelo
c
          if(abs(it).gt.imask) then
            ierr=15
            go to 90
          endif
          i=iand(it-iold+offset,mask)-offset
          iold=it
          is(x3+1)=i
          if(i.le.0) i=not(i)
          if(i.lt.64) then
c speed up for small differences
            y3=sig(i+2)
          else
            y3=8
            do while (y3.le.nb .and. i.ge.power(y3+1))
      	y3=y3+2	!was +1 (2/93, met)
            enddo
          endif
 7        y32=y3-y2
          if(y32.eq.0) go to 20	!cycle
          if(y12.eq.0) go to 18
          if(y12.lt.0) then
            if(y32.gt.0) then
c rise
      	x32=min(x3-x2,shift)
      	x20=min(x2,shift)
c low mid high --> mid high, low high
      	if(x20*y12+x32*y32.lt.0) then
      	  y2=y3
      	  y12=y1-y2
      	endif
            endif
c output if good enough
            if(y12*x2+shift.lt.0) go to 21
            y1=y2
            go to 18
          endif
c
          if(y32.lt.0) then
c drop
            if(y12*(x3-x2)+y12.gt.two) go to 21
c high mid low --> high low or high mid
            if(y12*(x3-x2)+y32.ge.0) go to 20
          else	!y32.gt.0
c dip
            if(y32.gt.y12) y32=y12
            y32=(x3-x2)*y32
            if(y32.gt.two) go to 21
            if(y32.gt.shift .and. x3.gt.256) go to 21
c high low mid --> high mid
c mid low high --> mid high
            if(y3-y1.lt.0) go to 19
          endif
 18       x2=x3
 19       y2=y3
          y12=y1-y2
 20     enddo
c
        if(y1.lt.0) then
          y1=y2
          x2=nn
        endif
 21     x2=min(x2,256)
        if((y1*x2+bit)/32.gt.maxp) then
          ierr=14
          go to 90
        endif
        if(y1.gt.16) then
          ierr=13
          go to 90
        endif
        ival=x2-1		!low order byte = count (max 255)
        ival=ior(ival,ishft(y1,8))	!high order byte = size
        ipak=pack(kptr)
        kpos=kpos-16
        itmp=0
        if(kpos.lt.0) itmp=ishft(ival,kpos+32)	!shift left
        if(kpos.ne.0) ival=ishft(ival,kpos)	!+=left, -=right
        if(kpos.gt.0) then
          ipak=ior(ipak,ival)
        else
          pack(kptr)=ior(ipak,ival)
          kptr=kptr+1
          ipak=itmp
          kpos=kpos+32
        endif
        pack(kptr)=ipak
        bit=bit+16
        if(y1.eq.0) go to 80
        kbits=y1
        do k=1,x2	!n=x2
          ival=is(k)
          ival=iand(ival,masks(kbits))	!drop sign extension
          kpos=kpos-kbits
          itmp=0
          if(kpos.lt.0) itmp=ishft(ival,kpos+32)	!shift left
          if(kpos.ne.0) ival=ishft(ival,kpos)	!+=left, -=right
          if(kpos.gt.0) then
            ipak=ior(ipak,ival)
          else
            pack(kptr)=ior(ipak,ival)
            kptr=kptr+1
            ipak=itmp
            kpos=kpos+32
          endif
          pack(kptr)=ipak
        enddo
        bit=bit+kbits*x2
 80     continue
        x0=x0+x2
c
c  extract byte pair (dmc)
        bytehi=bitems(1,x0)
        bytelo=bitems(2,x0)
c
       if(bytehi.lt.0) bytehi=bytehi+256
        if(bytelo.lt.0) bytelo=bytelo+256
        iold=256*bytehi+bytelo
c
      enddo
      ierr=0
c
      iwords=bit/32+1
      call packordr(pack,iwords,1)	! standardize byte order
c
 90   return
      end
c
      subroutine xpand2(pack,npack,bunpack,first,last,ierr)
c
c+expand compressed data, inverse of cmprs2
c
c   dmc 27 mar 1996 -- variant on cmprsi/xpandi, for portable binary
c     application.  see comments, cmprs2 suboutine
c
c pack	  is packed array
c bunpack  is unpacked array (byte stream, in byte pairs)
c first	  is index(1,n) of item to place in unpack(1)
c last	  is index of last item placed in unpack
c ierr	  is returned error code (0=ok)
c	  13 if  invalid data (internal nbits > 16)
c
c=======================================================================
c
      implicit integer(a-z)
c
      integer*1 bunpack(2,last)
c
      integer pack(npack)
c
      integer power(18)
c
      integer*4 masks(0:16)
c
      data power/
     >	  0,1,2,4,8,16,32,64,128,256,512,
     >    1024,2048,4096,8192,16384,32768,65536/
#ifdef __F90
      data masks/ Z"00000000",
     >	  Z"00000001",Z"00000003",Z"00000007",Z"0000000f",
     >	  Z"0000001f",Z"0000003f",Z"0000007f",Z"000000ff",
     >	  Z"000001ff",Z"000003ff",Z"000007ff",Z"00000fff",
     >	  Z"00001fff",Z"00003fff",Z"00007fff",Z"0000ffff"/
#else
      data masks/ '00000000'x,
     >	  '00000001'x,'00000003'x,'00000007'x,'0000000f'x,
     >	  '0000001f'x,'0000003f'x,'0000007f'x,'000000ff'x,
     >	  '000001ff'x,'000003ff'x,'000007ff'x,'00000fff'x,
     >	  '00001fff'x,'00003fff'x,'00007fff'x,'0000ffff'x/
#endif
c
c-----------------------------------------------------------------------
c
      lowmask=masks(8)
c
c  create local byte ordering for pack array
c
      call packordr(pack,npack,0)
c
      bit=0
      nb=16
c
      index=0
      ii=0
      value=0
      lower=0
      modulo=power(nb+2)
      upper=lower+modulo-1
      ierr=0
      do while(index.lt.last)
C				!get next 16 bits
        kptr=bit/32+1		!ptr into pack
        kpos=mod(bit,32)	!offset in word
        kpos=kpos+16
        if(kpos.lt.32) then
          ival=ishft(pack(kptr),-(32-kpos))	!shift right
          ival=iand(ival,masks(16))
        else
          kpos=kpos-32
          ival=iand(pack(kptr),masks(16-kpos))
          kptr=kptr+1
          if(kpos.gt.0) then
            ntx=ishft(pack(kptr),-(32-kpos))	!shift right
            ntx=iand(ntx,masks(kpos))
            ival=ior(ishft(ival,kpos),ntx)	!shift left
          endif
        endif
        bit=bit+16
        bits=ishft(ival,-8)		!high order byte = size
        bits=-bits
        n=iand(ival,lowmask)+1		!low order byte = count
        n=min(n,last-index)
        if(bits.gt.16) then
          ierr=13
          return
        endif
        if(bits.eq.0) then
          do k=1,n
          index=index+1
          if(index.ge.first) then
            ii=ii+1
            unpack=value
            bytehi=unpack/256
            bytelo=unpack-256*bytehi
            bunpack(1,ii)=bytehi
            bunpack(2,ii)=bytelo
          endif
          enddo
          go to 80
        endif
        kbits=abs(bits)
        do k=1,n
          kpos=kpos+kbits
          if(kpos.lt.32) then
            ival=ishft(pack(kptr),-(32-kpos))	!shift right
            ival=iand(ival,masks(kbits))
          else
            kpos=kpos-32
            ival=iand(pack(kptr),masks(kbits-kpos))
            kptr=kptr+1
            if(kpos.gt.0) then
      	ntx=ishft(pack(kptr),-(32-kpos))	!shift right
      	ntx=iand(ntx,masks(kpos))
      	ival=ior(ishft(ival,kpos),ntx)	!shift left
            endif
          endif
cx	    if(bits.lt.0) then	!bits is always negative
            imask=not(masks(kbits-1))		!extend sign bit
            if(iand(imask,ival).ne.0) ival=ior(ival,imask)
cx	    endif
          value=value+ival
          if(value.lt.lower) value=value+modulo
          if(value.gt.upper) value=value-modulo
          index=index+1
          if(index.ge.first) then
            ii=ii+1
            unpack=iand(value,masks(16))
            bytehi=unpack/256
            bytelo=unpack-256*bytehi
            bunpack(1,ii)=bytehi
            bunpack(2,ii)=bytelo
          endif
        enddo
        bit=bit+kbits*n
 80   continue
      enddo
c
 90   return
      end
c
C----------------------------------
C  new dmc 27 Mar 1996
C
C  standardize byte ordering of 4-wide byte stream
C
      subroutine packordr(pack,npack,idir)
C
      implicit NONE
C
      integer npack
      integer*4 pack(npack)  ! the data
      integer idir           ! =1 => standard, =0 => local ordering
C
      character*50 zdbg
C
      integer*4 jwork,ibyte,jbyte
      integer*1 bwork(4)
C
      integer ipack,ii,ishf
C
      integer*4, save :: masks(4)
      logical :: initflag = .FALSE.
      integer :: ibit,ibiti,indx
C
C  store in order:  hi byte ... lo byte
C
C  DMC -- get rid of non-portable code; use IBSET
C #ifdef __F90
C       data masks/Z"ff000000",Z"00ff0000",Z"0000ff00",Z"000000ff"/
C #else
C       data masks/'ff000000'x,'00ff0000'x,'0000ff00'x,'000000ff'x/
C #endif
C
      if(.not.initflag) then

         initflag = .TRUE.

         do ibit =0,31
            indx = 1 + ibit/8
            ibiti = 31-ibit

            masks(indx) = ibset(masks(indx),ibiti)
CCCC            write(6,*) ibit,ibiti,indx,masks
         enddo

      endif
C
C  loop through the data pack
C
      do ipack=1,npack
         if(idir.eq.1) then
C  standardize order  hi -> lo
            do ii=1,4
               ibyte=iand(pack(ipack),masks(ii))
               ishf=8*(ii-4)
               if(ishf.ne.0) then
                  jbyte=ishft(ibyte,ishf)
               else
                  jbyte=ibyte
               endif
               bwork(ii)=jbyte
            enddo
            call copyint4(bwork,pack(ipack))
         else
C  create locally desired ordering
            jwork=0
            do ii=1,4
               call byte4out(pack(ipack),ii,bwork(ii))
               ibyte=bwork(ii)
               if(ibyte.lt.0) ibyte=ibyte+256
               ishf=8*(4-ii)
               if(ishf.ne.0) then
                  jbyte=ishft(ibyte,ishf)
               else
                  jbyte=ibyte
               endif
               jwork=ior(jwork,jbyte)
            enddo
            pack(ipack)=jwork
         endif
      enddo
C
      return
      end
 
 
 
 
 
 
 
 
 
 
 
 
 
