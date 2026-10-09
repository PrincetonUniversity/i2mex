C---------------------------------------------------------------------
C  cftio
C
C
C  fortran shell for C-based binary i/o routines -- fixed length
C  512 byte blocks.
C
C---------------------------------------------------------------------
C
C  check that fortran "lun" is in range.
C
      subroutine cftlunck(ilun,sr,ier)
C
      character*(*) sr
C
      if((ilun.lt.0).or.(ilun.gt.139)) then
         ier=1
         write(6,1001) sr,ilun
 1001    format(' ?',a,' -- ilun=',i10,', argument out of range.')
      else
         ier=0
      endif
C
      return
      end
C---------------------------------------------------------------------
C
C call C routine to write a 512 block
C
      subroutine cftwrite(ilun,xstream,ier)
      use iso_c_binding, only: c_char
C
      character(kind=c_char) xstream(*)
C
      call cftlunck(ilun,'cftwrite',ier)
      if(ier.ne.0) return
C
      call cftcwr(ilun,xstream,ier)
C
      return
      end
C---------------------------------------------------------------------
C
C call C routine to read a 512 block
C
      subroutine cftread(ilun,xstream,ier)
      use iso_c_binding, only: c_char
C
      character(kind=c_char) xstream(*)
C
      call cftlunck(ilun,'cftread',ier)
      if(ier.ne.0) return
C
      call cftcrd(ilun,xstream,ier)
C
      return
      end
C---------------------------------------------------------------------
C
C call C routine to rewind a file
C
      subroutine cftrewnd(ilun,ier)
C
      call cftlunck(ilun,'cftread',ier)
      if(ier.ne.0) return
C
      call cftrwn(ilun)
C
      return
      end
C---------------------------------------------------------------------
C
C call C routine to close a file
C
      subroutine cftclose(ilun)
C
      call cftlunck(ilun,'cftread',ier)
      if(ier.ne.0) return
C
      call cftcls(ilun)
C
      return
      end
C
C---------------------------------------------------------------------
C
C call C routine to open a file
C
      subroutine cftopen(ilun,filnam,acod,iquiet,ier)
      use iso_c_binding, only: c_char, c_null_char
C
C  ilun must be btw 0 and 99
C
      integer ilun                      ! fortran like lun number
      character*(*) filnam              ! name of file
      character*(*) acod                ! access code, 'wb' for write
C                                       !              'rb' for read
      integer iquiet                    ! =1:  suppress common error messages
      integer ier                       ! completion code, 0=OK
C
C  local buffers
C
      character(kind=c_char) acbyts(20)           ! null terminated string
      character(kind=c_char) fnbyts(256)          ! null terminated string
C
C--------------------------------------------
C
      call cftlunck(ilun,'cftread',ier)
      if(ier.ne.0) return
C
 
      do ic=1,20
         acbyts(ic)=c_null_char
      enddo
C
      do ic=1,256
         fnbyts(ic)=c_null_char
      enddo
C
      ilfn=len(filnam)
      if(ilfn.gt.255) then
         write(6,1001) filnam(1:ilfn)
 1001    format(' ?cftopen:  filename too long:  '/1x,a)
         ier=2
         return
      else
         do ic=1,ilfn
            fnbyts(ic)=filnam(ic:ic)
         enddo
      endif
C
      ilac=len(acod)
      if(ilac.gt.19) then
         write(6,1002) acod(1:ilac)
         ier=3
 1002    format(' ?cftopen:  file access string too long:  ',a)
         return
      else
         do ic=1,ilac
            acbyts(ic)=acod(ic:ic)
         enddo
      endif
C
C  do it.
C
      call cftopn(ilun,fnbyts,acbyts,iquiet,ier)
C
      return
      end
