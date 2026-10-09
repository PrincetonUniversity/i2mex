subroutine trdatbuf_intcopy(str1,str2,i,i1,i2,intbuf,ivarsize,idest)

  character*(*), intent(in) :: str1,str2
  integer, intent(inout) :: i,i1,i2
  integer, intent(in) :: intbuf(*),ivarsize(*)
  integer, intent(inout) :: idest

  if(str1.eq.str2) then
     idest = intbuf(i1)
     i=i+1
     i1=i2+1
     i2=i2+ivarsize(i)
  endif

end subroutine trdatbuf_intcopy


subroutine trdatbuf_intcopy_a(str1,str2,i,i1,i2,intbuf,ivarsize,idest_a)

  character*(*), intent(in) :: str1,str2
  integer, intent(inout) :: i,i1,i2
  integer, intent(in) :: intbuf(*),ivarsize(*)
  integer, intent(inout) :: idest_a(*)

  if(str1.eq.str2) then
     idest_a(1:ivarsize(i)) = intbuf(i1:i2)
     i=i+1
     i1=i2+1
     i2=i2+ivarsize(i)
  endif

end subroutine trdatbuf_intcopy_a
