!.........................................................
!            SUBROUTINE TR_ADD_REPLACE
!.........................................................
 
subroutine tr_add_replace(zline,znew,zpat)
!
  use tr_getnl
  implicit NONE
!
  integer iline,jline,kline
  character*(*), intent(in) :: zline,znew,zpat
  character(120) :: sltext
! real, intent(in) :: rvar
!------------------------------
!
! print *,'In the routine tr_add_replace'
! print *,zline,znew
 
 
! If the pattern zline is in any string, replace it.
 
  do iline=1,mltext_nlines
   if (index(mltext(iline),zline).ne.0) then
     mltext(iline)=zline//'='//znew
   endif
  enddo
 
! If the pattern zline is not in any string, add a line.
! If zpat is not a null string, add the string containing
! zline after the line with the pattern zpat.
 
  jline=0
  do iline=1,mltext_nlines
   if (index(mltext(iline),zline).eq.0) then
     jline=jline+1
   endif
  enddo
 
  kline=0
  if (zpat.ne." ") then
   do iline=1,mltext_nlines
    if (index(mltext(iline),trim(zpat)).ne.0) then
     kline=iline
    endif
   enddo
  endif
 
  if((jline.eq.mltext_nlines).and.(kline.eq.0)) then
     mltext_nlines=mltext_nlines+1
     mltext(mltext_nlines)=zline//'='//znew
  elseif ((jline.eq.mltext_nlines).and.(kline.ne.0)) then
     sltext=zline//'='//znew
     do iline=mltext_nlines,kline+1,-1
      mltext(iline+1)=mltext(iline)
     enddo
     mltext(kline+1)=sltext
     mltext_nlines=mltext_nlines+1
  endif
 
  return
end subroutine tr_add_replace
