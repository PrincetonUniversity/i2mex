subroutine jetyr(nshot, cyear)
  !...   returns year ( string ) for JET shot number ; data from
  !      http://users.jet.efda.org/pages/codes-data/cpfweb/years/last100/index.htm
  !
  !     Author         J. E. Conboy ( James.Conboy@ccfe.ac.uk )
  !     Version        2.09
  !     Last Modified  04Apr2014
  !
  !     Modifications
  !
  !     04Apr2014      J. E. Conboy
  !                    Added 2012, 2013
  !
  !     22Feb2012      J. E. Conboy
  !                    Added ( 21 Jun ) 2011,   2012 ( after Iter wall shutdown )
  !     11Jan2011      J. E. Conboy
  !                    Still 79853..
  !     16Feb2010      J. E. Conboy
  !                    Last pulse 79853 on 23Oct09 ( followed by ITER wall shutdown )
  !     13Jan2008      J. E. Conboy
  !                    76329 was last pulse of 2008
  !     26Nov2007      J. E. Conboy
  !                    Last is 70750, no more plasma this year
  !     07/18/07       CLF: 
  !                    replaced end of yr 2005 = 64639 was 65653
  !     22Mar2007      R.V.Budny added last shot of 2006 
  !     24Apr2006      J. E. Conboy
  !                    Yr 2005 covers 63445 ( 05Apr2004 )
  !                                to 65653 ( 21Apr2006 )
  !     18Jul2007      CLF replaced by 64639 according to R. Budny
  !
  !                    All pulses labelled as 2005 are restart/commissioning 
  !                    or worse, & unsuitable for transp analysis 
  !
  !     18Jan2004      J. E. Conboy
  !                    Update last_shot(2003)
  !
  !     16Dec2003      J. E. Conboy
  !                    Rewritten
  !
  !     ================================================

  implicit none
  !
  integer, intent(in)             :: nshot          ! shot number
  character(len=*), intent(inout) :: cyear          ! year, as string '99' or '2002'
  integer                         :: n              ! local nshot
  character(len=4)                :: cy             ! local year string
  integer                         :: lyear          ! length of year
  integer                         :: iyear          ! year index
  integer, parameter              :: iyr0 = 1983    ! first known year
  integer, parameter              :: iyrl = 2020    ! last known year
  !
  integer, dimension(iyr0:iyrl)   :: last_shot= &
       (/     0,  3075,  6032, 11217, 14113, &
       18773, 21030, 23500, 27023, 27968, &
       27968, 33223, 35779, 39708, 43786, &
       46727, 49802, 52885, 54345, 57477, &
       62313, 63445, 64639, 69149, 70750, &
       76329, 79853, 79853, 81643, 83809, &
       85861, 87958, 89472, 92504, 92577, &
       93212, 96563, 99000  /) !last shot for 2020 not yet defined
  !
  n = max0( nshot, 1 )
  loop: do iyear=iyrl,iyr0,-1
    if( last_shot(iyear) < n )  exit loop
  end do loop
  !
  iyear = iyear + 1
  write(cy,'(i4.4)') iyear
  !
  lyear = len( cyear)
  if( lyear .gt. 3 )  then
    cyear( lyear-3:)  = cy
  else
    cyear( lyear-1:)  =  cy(3:4)
  end if
  !
  return
end subroutine jetyr
