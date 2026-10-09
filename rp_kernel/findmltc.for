      Subroutine FindMltC (Zline, Iptr, Mult, Cconst, Ier)
      use nltrdat_mod

C	In Character string Zline, starting at Zline(Iptr:Iptr) look
C	For constructs of the type "Mult*Cconst," or "Mult*Cconst "
C    	where Cconst has the form '        ' or '  '' '
C	Iptr points to "," after character string upon exit.
 
C	Arguments:
C	Input:
      Character*(*) Zline
      Integer       Iptr
 
C	Output:
C	Integer       Iptr - new pointer to after M*C  to next 1st '
      Integer       Mult
      Character*(*)  Cconst
      Integer       Ier

      Logical       Lfound
      Character*40  Cmult
      Character*1   Quote

      data Quote  /''''/     ! single '
C
C	------------------------------------------
      If (Zline(Iptr:Iptr) .Eq. '*') Then
          Ier = 1        ! * is first character - error.
          Go to 9999
      End If   ! *
 
      Lfound = .False.
      Ier = 0
      Mult = 1
      Istart = Iptr
 
      Do 100 I=Iptr, NrecLen
          If (Zline(I:I) .Eq. '*') Then
      	If (Lfound) Go To 9999   ! Error - 2 * in string
      	Cmult = Zline(Iptr:I-1)
      	Read (CMult, 9000, err=9999) Mult
 9000           Format(BN, I6)
      	Lfound = .True.    ! * found
      	Istart = I+1       ! Start after *
      	Go to 100
          End If   ! *
 
          If (Zline(I:I) .Eq. ','   .Or.
     1	        Zline(I:I) .Eq. Quote .Or.
     2          Zline(I:I) .Eq. ' '       ) Then    ! Termination of M*
 
      	  Go to 101
          End If   ! ,
 
  100   Continue
 
  101   Continue
 
C	Cmult = Zline(Istart:Iend)
C	Read (Cmult, 9001, Err=9999) Cconst
C9001   Format(BN, A)
 
C	Iptr = Iend+1
 
 
      Do 300 I=Istart,NrecLen
          If (Zline(I:I) .Eq. Quote ) Then
C	      .First Quote found.
            I1 = I+1
            Do 200 II=I1,Nreclen
 
              If (Zline(II:II) .Eq. Quote) Then
C		  .Closing quote (#2) found.
      	  I2 = II-1
      	  Iptr = II+1   ! Point to trailing ","
 
      	  If (Zline(Iptr:iptr) .Ne. ','  .And.
     1                ZLine(Iptr:Iptr) .Ne. ' '       ) Then
      	      Ier = 9     ! No trailing blank or ","
      	      Go To 9999
                End If    ! ","
 
      	  Go To 301
      	End If   ! '
  200         Continue
            Ier = 1   ! No second quote
            Go to 9999
          End If      ! '1
 
  300     Continue
        Ier = 1    ! No first quote
        Go to 9999
 
  301   Continue
      CConst = Zline(I1:I2)
 
      Return
 
 
 9999 Continue   ! error exit
      Ier = 1
      Return
      End
