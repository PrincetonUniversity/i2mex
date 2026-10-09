      Subroutine FindMltR (Zline, Iptr, Mult, Rconst, Ier)
      use nltrdat_mod

C	In Character string Zline, starting at Zline(Iptr:Iptr) look
C	For constructs of the type "Mult*Rconst," or "Mult*Rconst "
 
C	Arguments:
C	Input:
      Character*(*) Zline
      Integer       Iptr
 
C	Output:
C	Integer       Iptr - new pointer to after M*R,
      Integer       Mult
      Real          Rconst
      Integer       Ier
 
      Logical       Lfound
      Character*40  Cmult
 
C	------------------------------------------
 
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
 
          If (Zline(I:I) .Eq. ',' .Or.
     1          Zline(I:I) .Eq. ' '      ) Then    ! Termination of M*R or R
      	  Iend = I-1
 
      	  Go to 101
          End If   ! ,
 
  100   Continue
      Iend = NrecLen
 
  101   Continue
 
      Cmult = Zline(Istart:Iend)
      Read (Cmult, 9001, Err=9999) Rconst
 9001   Format(BN, E15.0)
 
      Iptr = Iend+1
      Return
 
 
 9999 Continue   ! error exit
      Ier = 1
      Return
      End
