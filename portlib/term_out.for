C***********************************************************************
      SUBROUTINE TERM_OUT(BUFFER,NBYTES)
C
      use iso_c_binding, only: c_char
      implicit none
C
C
C...This routine was written by Harry H. Towner of the Princeton
C...Plasma Physics Lab. This routine is similar to John Coonrod's
C...routine JC_LINE. This routine will output NBYTES of array BUFFER
C...to the user's terminal
C
C...Parameters:
C*..BUFFER	- Array of NBYTES to be sent to ther user's terminal.
C*..NBYTES	- The number of bytes to be sent to the terminal
C
C***********************************************************************
C
      integer nbytes
      character(kind=c_char) BUFFER(NBYTES)
C
      character(kind=c_char) :: itemp
      INTEGER ISTAT,I,ZPUTC
C
      DO I=1,NBYTES
        ITEMP=BUFFER(I)
        ISTAT=ZPUTC(ITEMP)
      END DO
C
      RETURN
      END
C-----------------------
      subroutine term_str_out(str)
      use iso_c_binding, only: c_char
      implicit none
C
      character*(*) str
C
      integer NBSIZ
      parameter (NBSIZ=128)
C
      character(kind=c_char) BUFFER(NBSIZ)
C
      integer ilen,ic0,ice,inb,i,ic
C
      ilen=len_trim(str)
C
      do ic0=1,ilen,NBSIZ
         ice=min(ilen,(ic0+NBSIZ-1))
         inb=ice-ic0+1
         do i=1,inb
            ic=ic0+i-1
            buffer(i)=str(ic:ic)
         enddo
C
         call term_out(buffer,inb)
      enddo
C
      return
      end
