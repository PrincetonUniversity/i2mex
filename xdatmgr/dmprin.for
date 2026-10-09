C******************** START FILE DMPRIN.FOR ; GROUP PLDMGR ******************
C-----------------------------------------------------------
C  DMPRIN
C    DEBUG... PRINTOUT CONTENTS OF DATA AREA
C
C    write report to unit 6
C
      SUBROUTINE DMPRIN(str,jmax)
 
      implicit NONE

      character*(*), intent(in) :: str  ! label string: caller
      integer, intent(in) :: jmax       ! max no. of entries to display

C  Updated:
C
C  05/01/11 dmc break into subroutines; allow output to files
C     DMPRIN -> DMPRINU, DMPRIN_FILE, DMPRIN_GF
C
C  04/17/09 dmc added str argument
C   1/20/94 tbt Added counter J.
C  09/28/93 TBT Added Lacc & Lnext to output.
C
      call dmprinu(6,str,jmax)
C
      END
C----------------------------------------------------
      SUBROUTINE DMPRINU(ilun,str,jmax)
C
C    write report to unit (ilun)
C
      use datmgr_mod
      implicit NONE
C
      integer, intent(in) :: ilun       ! Fortran unit to write
      character*(*), intent(in) :: str  ! label string: caller
      integer, intent(in) :: jmax       ! max no. of entries to display

C----------------------------------------------
      Integer J,JL
      character*2 suffix
C----------------------------------------------

      write(ilun,*) ' ----------> DMPRIN call from: ',trim(str)

      J =0
      JL=1
 10   CONTINUE
      J =J+1
      If (J .GT. jmax) Go to 100    ! Prevent infinite loop. tbt
 
      suffix=' '
      if(lnext(jl).gt.0) then
         if(lprev(lnext(jl)).ne.jl) suffix(1:1)='?'
      endif

      if(lprev(jl).gt.0) then
         if(lnext(lprev(jl)).ne.jl) suffix(2:2)='!'
      endif

      write(ilun,1001) J, JL,
     1             DMGLBL(JL),LOCD(JL),NWDS(JL),MPRIO(JL),
     1             LACC(JL), LNEXT(JL), LPREV(JL),suffix
 
 1001 FORMAT(' ',I3, I3, ' "',A,'"'/
     1   '        LOC=',I9,' SIZE=',I9,' PRIO=',I2,
     1   ' LACC=', I5, ' Lnext=', I4, ' Lprev=', I4,1x,A)
 
      IF(DMGLBL(JL).EQ.'%FINI') GO TO 100
      JL=LNEXT(JL)
      GO TO 10
C
 100  CONTINUE
C
      RETURN
      END

C----------------------------------------------------
      SUBROUTINE DMPRIN_FILE(filename,str,jmax)

C  write report to file, as specified...

      character*(*), intent(in) :: filename
      character*(*), intent(in) :: str  ! label string: caller
      integer, intent(in) :: jmax       ! max no. of entries to display

      !----------
      integer :: io
      !----------

      call find_io_unit(io)

      open(unit=io,file=filename,status='unknown')

      call dmprinu(io,str,jmax)

      close(unit=io)

      END

C----------------------------------------------------
      SUBROUTINE DMPRIN_GF(ict,str,jmax)

C  generate filename from integer, write report to file

      integer, intent(in) :: ict        ! input integer: 1,2,3... < 10000
      character*(*), intent(in) :: str  ! label string: caller
      integer, intent(in) :: jmax       ! max no. of entries to display

C------------------------------
      character*30 :: filename
C------------------------------

      filename = ' '
      write(filename,'("dmprin",i4.4,".dat")') ict

      call dmprin_file(filename,str,jmax)

      END
