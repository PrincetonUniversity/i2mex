      subroutine tfsclab(ilun,zlbl,inum)

      use cplotr_mod

      integer, intent(in) :: ilun    ! unit to write
      character*(*), intent(in) :: zlbl  ! passed label
      integer, intent(in) :: inum    ! scalar function number

      !  write out full label of a scalar function

      if((inum.lt.1).or.(inum.gt.nft)) then
         write(ilun,*) trim(zlbl)//':  index inum=',inum
         write(ilun,*) '    out of range: 1 to nft=',nft
      else
         write(ilun,*) trim(zlbl)//': '//abt(inum)
         write(ilun,*) '   ',labelt(inum),'   ',unitst(inum)
      endif

      return
      end
