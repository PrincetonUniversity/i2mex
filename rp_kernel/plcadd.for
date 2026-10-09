      Subroutine PlcAdd (A, B)
 
C	.Created 11/6/92 by Tbt. Called from PlcFnxct in RPLOT to get
C	.around a subscript problem which occurs in the handler when
C	.an arithmetic exception happens.
 
C	.Add two numbers.
 
      Double Precision A,B
 
      A = A + B
	
      Return
      End
