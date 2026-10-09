!******************** START FILE FPOLAR.FOR ; GROUP TKBNDRY ******************
!.....................................................................
!  REV DMC -- REMOVED TRANSP COMMON
!   SO I CAN USE THIS IN RPLOT
!
!  REV BB / MAY 94 -- USE ATAN2 INSTEAD OF ATAN TO EVALUATE FPOLAR
!
!   output in range [0,twopi] instead of range [-pi,pi]
!   atan2 and fpolar are the same in the upper half plane; fpolar
!   is twopi greater in the lower half plane.
!
      REAL FUNCTION FPOLAR ( ZR77, ZZ77 )
!
!	CALC. THE POLAR ANGLE DEFINED BY TAN(THETA)=ZZ77/ZR77
!
      DATA TWOPI/6.2831853071795862E+00/
      DATA PI/3.1415926535897931E+00/
!
!--------------------------------------------
!
      if (zz77 .ge. 0.0) then
!
!  upper half plane  --  including z=0 line
!
       if (ZR77 .NE. 0.) THEN
!
          FPOLAR = ATAN2 ( ZZ77 , ZR77 )
!
       else
!
!  on the vertical axis, at or above the z=0 line
!
        FPOLAR = PI*0.5
!
       end if
!
      else
!
!  lower half plane  --  excluding z=0 line
!
          FPOLAR = ATAN2 ( ZZ77 , ZR77 ) + TWOPI
!
      end if
!
      return
      end
!******************** end FILE FPOLAR.FOR ; GROUP TKBNDRY ******************
