subroutine mastyr( nshot, cyear )
  !...   returns year ( string ) for MAST shot number ; data from
  !      Carol Brickley
  !      (or see http://fusweb1.fusion.ccfe.ac.uk/~mast_sl/SLlog/sllog_index.html,
  !       on CCFE intranet )
  !
  !     Author         J. E. Conboy ( James.Conboy@jet.uk )
  !     Version        2.03
  !     Last Modified  11Jan2011
  !
  !     Modifications
  !
  !     23Apr2014      jim.conboy@ccfe.ac.uk
  !                    Added last shot for 2013 ;  followed by MAST-U upgrade
  !
  !     22Feb2013      C. Ludescher-Furth
  !                    The 2012 range was 27963-28410 inclusive
  !     22Feb2012      J. E. Conboy
  !                    The 2011 range was 25402-27962 inclusive
  !     11Jan2011      J. E. Conboy
  !                    25401 is last shot of 2010
  !     16Feb2010      J. E. Conboy
  !                    23944 is last shot of 2009
  !     13Jan2009      J. E. Conboy
  !                    21089 is last shot of 2008
  !     16Dec2003      J. E. Conboy
  !                    Cribbed from jetyr
  !
  !     ================================================

  implicit none
  !
  integer, intent(in)             :: nshot         ! shot number
  character(len=*), intent(inout) :: cyear         ! year, as string - '99' or '2002'
  integer                         :: n             ! local nshot
  character(len=4)                :: cy            ! local year string
  integer                         :: lyear         ! length of year
  integer                         :: iyear         ! year index
  integer, parameter              :: iyr0 = 1999   ! first known year
  integer, parameter              :: iyrl = 2013   ! last known year
  !
  integer,dimension(iyr0:iyrl)    :: last_shot= &
       (/     0, &
       3416,  5009,  7951,  9989, &
       11725, 14706, 16660, 19026, 21089, &
       23944, 25401, 27962, 28410, 30473 /)
  !
  n = max0( nshot, 1 )
  loop:  do iyear=iyrl,iyr0,-1
    if( last_shot(iyear) < n ) exit loop
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
end subroutine mastyr
