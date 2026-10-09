      Subroutine FindName ( Varname, Lenv, Ind, Ier)
      use nltrdat_mod

C	Find Varname in the list of Namelist names.
C   	Return its index in that list in Ind.
C	Start search at index Ind.
C	Return the length of the name found in Lenv
 
C	Arguments:
C	Input:
      Character*(*) Varname
      Integer       Ind
 
C	Output:
C	Integer       Ind
      Integer       Ier
 
      If (Ind .Le. 0  .Or. Ind .Gt. NNames) Go to 901  ! Error
      Ier = 0
 
      Do 100 I=Ind,NNames
          If (Varname .Eq. Names(I))  goto 102
  100   Continue
 
C	.Name not found - error return
      Go to 901
 
 
  102 Continue   ! Name found
      Ind = I
 
      Lenvar = Len(varname)
      Do 200 I=Lenvar,1,-1
          If (Varname(I:I) .Ne. ' ') Go to 201
  200   Continue
 
C	.Err - name all blank  !! impossible
      Go to 901
 
  201   Continue
      Lenv = I   ! Last non-blank lettter
 
      Return
 
 
  901   Continue
      ier = 1
      Ind = 0
      Lenv = 0
      Return
 
      End
