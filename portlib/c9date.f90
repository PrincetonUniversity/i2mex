subroutine c9date(zdate)
  implicit none
  character, intent(out) :: zdate*9
  character :: zbuf*24
!
!  get a 9 character ddmmmyyyy
!
  call fdate(zbuf)
  if(zbuf(9:9).eq.' ')zbuf(9:9)=char(48)
  zdate=zbuf(9:10)//zbuf(5:7)//zbuf(21:24)
  return
end subroutine c9date
