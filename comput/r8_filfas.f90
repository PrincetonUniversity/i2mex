!******************** START FILE FILFAS.FOR ; GROUP FILFAS ******************
!---------------------------------------------------------------
!  FILFAS   D. Mc Cune  pppl 4 dec 1981
!  A FAST MOVING TRIANGULAR WEIGHTED AVERAGE SMOOTHING ROUTINE
!
subroutine R8_FILFAS(F,FSM,N,NWA)
  use iso_c_binding, only: fp => c_double
  implicit none
  integer :: n,nwa,inwa,j,k,iw,ipt1,ipt2
  real(fp), dimension(n) :: fsm,f
  real(fp) :: fac,sum
  !
  !  INPUT:
  !   F(N) DATA TO BE SMOOTHED
  !   NWA-- NUMBER OF POINTS ON EITHER SIDE OF F(J) TO AVERAGE
  !         WITH TRIANGULAR WEIGHTING TO COMPUTE FSM(J)
  !
  !  OUTPUT:
  !   FSM(N): SMOOTHED DATA
  !
  !
  !  ALGORITHM:
  !  FSM(J) IS A WEIGHTED AVERAGE OF F(J-NWA),F(J-NWA+1),...,
  !                                  F(J),...,F(J+NWA-1),F(J+NWA)
  !  THE WEIGHTING FACTOR FOR
  !  F(J+K) IS W(K)=(NWA+1-IABS(K))/(NWA+1)**2
  !   NOTE THAT SUM[K= -NWA TO +NWA]W(K) = 1.0 EXACTLY
  !
  !  W(K) IS A STEP FUNCTION APPROXIMATION TO A TRIANGULAR WEIGHTING
  !  FUNCTION
  !
  !  END CONDITIONS:  if THE INDEX (J+K) IS .LT. 1 IT IS
  !   REPLACED BY 1; if (J+K) .GT. N IT IS REPLACED BY N.
  !  THUS THE END POINTS ARE MORE HEAVILY WEIGHTED IN THE REGIONS NEAR THE
  !  ENDPOINTS OF THE DATA BEING SMOOTHED
  !
  !  EXECUTABLE CODE
  !
  if(NWA.LE.0) GO TO 200
  INWA=NWA+1
  FAC=1.0D0/(INWA**2)
  !  START LOOP OVER DATA POINTS
  do J=1,N
    !  INITIATE SUM FOR POINT J
    SUM=F(J)*INWA
    !  LOOP OVER POINTS TO AVERAGE
    do K=1,NWA
      IW=INWA-K
      IPT1=max(1,(J-K))
      IPT2=min(N,(J+K))
      SUM=SUM+IW*(F(IPT1)+F(IPT2))
    end do
    FSM(J)=SUM*FAC
  end do
  return
  !  NWA < 0: JUST COPY DATA
200 continue
  do J=1,N
    FSM(J)=F(J)
  end do
  return
end subroutine R8_FILFAS
