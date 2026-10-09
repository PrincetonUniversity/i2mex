C******************** START FILE CAPTL.FOR ; GROUP SPLITN ******************
c-----------------------------------------------------
c  capitolize all characters in input line.   RPLOT version. TBT 3/93
c
      subroutine captl(line)
c
      character*(*) line
c
      character*52 alphab
      data alphab/
     >  'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'/
c
      ilen=len(line)
      if(ilen.le.0) return
 
      do 20 i=1,ilen
         do 10 j=1,26
           if(line(i:i).eq.alphab(j:j)) Then
              line(i:i)=alphab(j+26:j+26)
              Go to 20
           End if    ! line
 10      continue
 20   Continue
 
      return
      end
C******************** END FILE CAPITL.FOR ; GROUP SPLITN ******************
