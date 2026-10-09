subroutine splitn_printwarn(msg,nlist,iwid,list,klist)

  !  if list is not empty, print message and display affected items

  implicit NONE
  character*(*), intent(in) :: msg   ! message to print (no trailing blanks)
  integer, intent(in) :: nlist       ! length of list
  integer, intent(in) :: iwid        ! width of list elements
  character*(iwid), intent(in) :: list(nlist)  ! the list
  integer, intent(in) :: klist(nlist)   ! block affiliation of list

  !---------------------------------------
  character*78 buf
  integer i,j,jb,toggle,maxblock
  character*6 :: bencod
  !---------------------------------------

  if(nlist.eq.0) return

  write(6,*) ' ' 
  write(6,*) msg

  maxblock=0
  toggle=0

  do i=1,nlist
     if(toggle.eq.0) then
        if(klist(i).gt.0) then
           call mk_bencod
           buf='  '//trim(list(i))//'['//bencod(jb:6)//']'
        else
           buf='  '//list(i)
        endif
     else
        if(klist(i).gt.0) then
           call mk_bencod
           buf(40:)=trim(list(i))//'['//bencod(jb:6)//']'
        else
           buf(40:)=list(i)
        endif
        write(6,*) buf
     endif
     toggle=1-toggle   ! 0,1,0,1,...
  enddo

  if(toggle.eq.1) write(6,*) buf
  if(maxblock.gt.0) then
     write(6,*) &
          '  "[ ]" indicates warning pertaining to indicated update block.'
  endif

CONTAINS

  subroutine mk_bencod

    maxblock = max(maxblock,klist(i))
    bencod=' '
    write(bencod,'(i6)') klist(i)
    jb=1
    do j=5,1,-1
       if(bencod(j:j).eq.' ') then
          jb=j+1
          exit
       endif
    enddo

  end subroutine mk_bencod

end subroutine splitn_printwarn
