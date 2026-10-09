subroutine splitn_chkdp(zlin1,zlin2)

  use splitn_module
  implicit NONE

  character*(*), intent(inout) :: zlin1  ! namelist file input line (orig.)
  character*(*), intent(inout) :: zlin2  ! namelist file input line (uppercase)

  !  check that all non-null floating point value fields include a decimal
  !  point (needed for reliable portable reads).  If a value field lacks a
  !  decimal point, then, insert it in both zlin1 and zlin2, and increment
  !  the module's affected line character indices:  val2  lcmt  len1.

  !---------------------------
  integer iscan
  integer iw1,iw2,ir1,ir2,idelim,irepeat,iextra
  integer idp,insrt,inam2
  character*150 ztmp
  !---------------------------

  iscan = val1 - 1

  do
     call splitn_nxgrp(zlin2,iscan,iw1,iw2,idelim,ir1,ir2,irepeat,iextra)
     if(iw2.ge.iw1) then
        !  look for decimal pt.
        idp=index(zlin2(iw1:iw2),'.')

        if(idp.eq.0) then
           !  ...not found: insert & add to warning list
           !     decimal pt to be inserted before exponent field, or,
           !     at the end if there is not exponent field.

           inam2=index(zlin2(nam1:nam2),'(')
           if(inam2.eq.0) then
              inam2=nam2          ! point to end of name
           else
              inam2=inam2+nam1-2  ! point before left parenthesis
           endif
           call splitn_addwarn(zlin2(nam1:inam2),lackdp,nlackdp, &
                len(lackdp(1)),klackdp,maxwarn)

           insrt=index(zlin2(iw1:iw2),'E')
           if(insrt.eq.0) insrt=index(zlin2(iw1:iw2),'D')
           if(insrt.eq.0) then
              insrt=iw2+1
           else
              insrt=insrt+iw1-1
           endif

           ztmp=zlin1(1:insrt-1)//'.'
           ztmp(insrt+1:)=zlin1(insrt:)
           zlin1=ztmp

           ztmp=zlin2(1:insrt-1)//'.'
           ztmp(insrt+1:)=zlin2(insrt:)
           zlin2=ztmp

           iw2=iw2+1
           iscan=iscan+1
           idelim=idelim+1

           len1=len1+1
           if(lcmt.gt.0) lcmt=lcmt+1
           val2=val2+1

        endif

     endif

     if(idelim.ge.val2) then
        exit
     else
        iscan=idelim
     endif
  enddo

end subroutine splitn_chkdp
