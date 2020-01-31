
subroutine gaussPointsAndWeights(nt1, t, ngauss, tgauss, wgauss, ier)
  
  ! Given a grid t and 1 >= ngauss <=10, return the Gauss 
  ! quadrature points tgauss and weights wgauss. *Note* tgauss and wgauss
  ! must be preallocatted and must have size (nt1-1)*ngauss.
  !
  ! pletzer Fri Sep  8 15:35:59 EDT 2000

  implicit none
  integer, parameter :: r8 = selected_real_kind(12,100)

  integer, intent(in) :: nt1 ! # of grid t-points
  real(r8), intent(in) :: t(*)  ! grid
  integer, intent(in) :: ngauss ! number of Gauss points per interval
  real(r8), intent(out) :: tgauss(*) ! the output array of t-Gauss points
  real(r8), intent(out) :: wgauss(*) ! the output array of Gauss weights
  integer, intent(out) :: ier ! error flag: 0=ok, 
                              !             1=ngauss too big, 
                              !             2=ngauss too small
                              !             3=wrong size of tgauss
                              !             4=wrong size of wgauss

  integer, parameter :: ngmax = 10
  real(r8) :: dt
  real(r8), dimension(ngmax) :: point, weigh
  integer i,j


  ier = 0
!!$  if(size(tgauss) < ngauss*(nt1-1)) then
!!$     ier = 3
!!$     print *,'--gaussPointsAndWeights-- wrong size of tgauss ', size(tgauss)
!!$     return
!!$  endif
!!$  if(size(wgauss) < ngauss*(nt1-1)) then
!!$     ier = 4
!!$     print *,'--gaussPointsAndWeights--wrong size of wgauss ', size(wgauss)
!!$     return
!!$  endif

 SELECT CASE(ngauss)


    case (1)  ! 1 Gauss' points

      point(1) = 0._r8;
 
      weigh(1) = 2._r8;
 
 
    case (2)  ! 2 Gauss' points

      point(1) = -0.5773502691896257_r8;
      point(2) = 0.5773502691896257_r8;
 
      weigh(1) = 1._r8;
      weigh(2) = 1._r8;
 
 
    case (3)  ! 3 Gauss' points

      point(1) = -0.7745966692414834_r8;
      point(2) = 0._r8;
      point(3) = 0.7745966692414834_r8;
 
      weigh(1) = 0.5555555555555553_r8;
      weigh(2) = 0.888888888888889_r8;
      weigh(3) = 0.5555555555555553_r8;
 
 
    case (4)  ! 4 Gauss' points

      point(1) = -0.861136311594053_r8;
      point(2) = -0.3399810435848562_r8;
      point(3) = 0.3399810435848562_r8;
      point(4) = 0.861136311594053_r8;
 
      weigh(1) = 0.3478548451374539_r8;
      weigh(2) = 0.6521451548625462_r8;
      weigh(3) = 0.6521451548625462_r8;
      weigh(4) = 0.3478548451374539_r8;
 
    case (5)  ! 5 Gauss' points

      point(1) = -0.906179845938664_r8;
      point(2) = -0.5384693101056829_r8;
      point(3) = 0._r8;
      point(4) = 0.5384693101056829_r8;
      point(5) = 0.906179845938664_r8;
 
      weigh(1) = 0.2369268850561887_r8;
      weigh(2) = 0.4786286704993665_r8;
      weigh(3) = 0.5688888888888889_r8;
      weigh(4) = 0.4786286704993665_r8;
      weigh(5) = 0.2369268850561887_r8;
 
 
    case (6)  ! 6 Gauss' points

      point(1) = -0.932469514203152_r8;
      point(2) = -0.6612093864662646_r8;
      point(3) = -0.2386191860831968_r8;
      point(4) = 0.2386191860831968_r8;
      point(5) = 0.6612093864662646_r8;
      point(6) = 0.932469514203152_r8;
 
      weigh(1) = 0.1713244923791709_r8;
      weigh(2) = 0.3607615730481379_r8;
      weigh(3) = 0.4679139345726913_r8;
      weigh(4) = 0.4679139345726913_r8;
      weigh(5) = 0.3607615730481379_r8;
      weigh(6) = 0.1713244923791709_r8;
 
 
    case (7)  ! 7 Gauss' points

      point(1) = -0.949107912342759_r8;
      point(2) = -0.7415311855993937_r8;
      point(3) = -0.4058451513773972_r8;
      point(4) = 0._r8;
      point(5) = 0.4058451513773971_r8;
      point(6) = 0.7415311855993945_r8;
      point(7) = 0.949107912342759_r8;
      weigh(1) = 0.129484966168868_r8;
      weigh(2) = 0.2797053914892783_r8;
      weigh(3) = 0.3818300505051186_r8;
      weigh(4) = 0.4179591836734694_r8;
      weigh(5) = 0.3818300505051188_r8;
      weigh(6) = 0.279705391489276_r8;
      weigh(7) = 0.1294849661688697_r8;
 
 
    case (8)  ! 8 Gauss' points

      point(1) = -0.960289856497537_r8;
      point(2) = -0.7966664774136262_r8;
      point(3) = -0.5255324099163289_r8;
      point(4) = -0.1834346424956498_r8;
      point(5) = 0.1834346424956498_r8;
      point(6) = 0.5255324099163289_r8;
      point(7) = 0.7966664774136262_r8;
      point(8) = 0.960289856497537_r8;
 
      weigh(1) = 0.1012285362903738_r8;
      weigh(2) = 0.2223810344533786_r8;
      weigh(3) = 0.3137066458778874_r8;
      weigh(4) = 0.3626837833783619_r8;
      weigh(5) = 0.3626837833783619_r8;
      weigh(6) = 0.3137066458778874_r8;
      weigh(7) = 0.2223810344533786_r8;
      weigh(8) = 0.1012285362903738_r8;
 
 
    case (9)  ! 9 Gauss' points

      point(1) = -0.968160239507626_r8;
      point(2) = -0.836031107326637_r8;
      point(3) = -0.6133714327005903_r8;
      point(4) = -0.3242534234038088_r8;
      point(5) = 0._r8;
      point(6) = 0.3242534234038088_r8;
      point(7) = 0.6133714327005908_r8;
      point(8) = 0.836031107326635_r8;
      point(9) = 0.968160239507627_r8;
 
      weigh(1) = 0.0812743883615759_r8;
      weigh(2) = 0.1806481606948543_r8;
      weigh(3) = 0.2606106964029356_r8;
      weigh(4) = 0.3123470770400029_r8;
      weigh(5) = 0.3302393550012597_r8;
      weigh(6) = 0.3123470770400025_r8;
      weigh(7) = 0.2606106964029353_r8;
      weigh(8) = 0.1806481606948577_r8;
      weigh(9) = 0.0812743883615721_r8;
 
 
    case (10)  ! 10 Gauss' points

      point(1) = -0.973906528517172_r8;
      point(2) = -0.865063366688984_r8;
      point(3) = -0.6794095682990246_r8;
      point(4) = -0.433395394129247_r8;
      point(5) = -0.1488743389816312_r8;
      point(6) = 0.1488743389816312_r8;
      point(7) = 0.433395394129247_r8;
      point(8) = 0.6794095682990246_r8;
      point(9) = 0.865063366688984_r8;
      point(10) = 0.973906528517172_r8;
 
      weigh(1) = 0.06667134430868681_r8;
      weigh(2) = 0.149451349150573_r8;
      weigh(3) = 0.2190863625159832_r8;
      weigh(4) = 0.2692667193099968_r8;
      weigh(5) = 0.2955242247147529_r8;
      weigh(6) = 0.2955242247147529_r8;
      weigh(7) = 0.2692667193099968_r8;
      weigh(8) = 0.2190863625159832_r8;
      weigh(9) = 0.149451349150573_r8;
      weigh(10) = 0.06667134430868681_r8;
    
   case (ngmax+1:)
      
      ier = 1
      print *,'--gaussPointsAndWeights--ngauss=',ngauss,' too big'
      return

   case (:0)

      ier = 2
      print *,'--gaussPointsAndWeights--ngauss=',ngauss,' too small'
      return
 
 END SELECT
  

 do i = 1, nt1 -1
    dt = t(i+1) - t(i)
    do j = 1, ngauss
       tgauss((i-1)*ngauss + j) = t(i) + dt*(point(j) + 1.0_r8)/2.0_r8
       wgauss((i-1)*ngauss + j) = dt*weigh(j)/2.0_r8
    enddo
 enddo


end subroutine 


 
