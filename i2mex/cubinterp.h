!*COMDECK CUCCCC
! ----------------------------------------------------------------------
! --     STATEMENT FUNCTION FOR CUBIC INTERPOLATION                   --
! --                         23.04.88            AR        CRPP       --
! --                                                                  --
! -- CUBIC INTERPOLATION OF A FUNCTION F(X)                           --
! -- THE EIGHT ARGUMENTS A1,A2,A3,A4,B1,B2,B3,B4 ARE DEFINED BY:      --
! -- F(B1) = A1 , F(B2) = A2 , F(B3) = A3 , F(B4) = A4                --
! ----------------------------------------------------------------------
!
 REAL*8 FA3, FA2, FA1, FA0
 REAL*8 FCCCC0
 REAL*8 A1,A2,A3,A4,B1,B2,B3,B4,PX
!
         FA3(A1,A2,A3,A4,B1,B2,B3,B4) = &
     &        (A1-A2) / ((B1-B2)*(B2-B4)*(B2-B3)) + &
     &        (A1-A3) / ((B4-B3)*(B3-B1)*(B3-B2)) + &
     &        (A1-A4) / ((B1-B4)*(B2-B4)*(B3-B4))
         FA2(A1,A2,A3,A4,B1,B2,B3,B4) = &
     &        (A1-A2) / ((B2-B1)*(B3-B2)) + &
     &        (A3-A1) / ((B3-B1)*(B3-B2)) - &
     &        (B1+B2+B3) * FA3(A1,A2,A3,A4,B1,B2,B3,B4)
         FA1(A1,A2,A3,A4,B1,B2,B3,B4) = &
     &        (A1-A2) / (B1-B2) - &
     &        (B1+B2) * FA2(A1,A2,A3,A4,B1,B2,B3,B4) - &
     &        (B1*B1+B1*B2+B2*B2) * FA3(A1,A2,A3,A4,B1,B2,B3,B4)
         FA0(A1,A2,A3,A4,B1,B2,B3,B4) = &
     &        A1 - &
     &        B1 * (FA1(A1,A2,A3,A4,B1,B2,B3,B4) + &
     &              B1 * (FA2(A1,A2,A3,A4,B1,B2,B3,B4) + &
     &                    B1 * FA3(A1,A2,A3,A4,B1,B2,B3,B4)))
! ----------------------------------------------------------------------
! -- FCCCC0 GIVES THE VALUE OF THE FUNCTION AT POINT PX:              --
! -- FCCCC0(......,PX) = F(PX)                                        --
! ----------------------------------------------------------------------
        FCCCC0(A1,A2,A3,A4,B1,B2,B3,B4,PX) = &
     &              FA0(A1,A2,A3,A4,B1,B2,B3,B4) + &
     &              PX * (FA1(A1,A2,A3,A4,B1,B2,B3,B4) + &
     &                    PX * (FA2(A1,A2,A3,A4,B1,B2,B3,B4) + &
     &                          PX * FA3(A1,A2,A3,A4,B1,B2,B3,B4)))
