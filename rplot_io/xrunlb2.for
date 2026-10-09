      subroutine xrunlb2(runlb2,part1,ilp1,part2,ilp2,part3,ilp3)
c
c  break legacy runlb2 label "tok.yy runid (MDS+)" into its 3
c  constituent parts -- use blank delimiter
c
      character*(*) runlb2              ! input run label
c
c  output:
c
      character*(*) part1               ! 1st word
      integer ilp1                      ! non-blank length
c
      character*(*) part2               ! 2nd word
      integer ilp2                      ! non-blank length
c
      character*(*) part3               ! 3rd word
      integer ilp3                      ! non-blank length
c
c  typically:
c    part1 = tok.yy, e.g. "TFTR.96"
c    part2 = runid, e.g. 12345A06
c    part3 = blank or "(MDS+)"
c
c------------------------------------------
c
      integer str_length
c
c------------------------------------------
c
      part1=' '
      part2=' '
      part3=' '
c
      ilp1=1
      ilp2=1
      ilp3=1
c
c  find first non-blank in runlb2
c
      ilnb=str_length(runlb2)
c
      do ic=1,ilnb
         if(runlb2(ic:ic).ne.' ') go to 10
      enddo
c
      return                            ! all blanks...
c
 10   continue
      icp1a=ic
c
      do ic=icp1a+1,ilnb
         if(runlb2(ic:ic).eq.' ') go to 20
      enddo
      ic=ilnb+1
c
 20   continue
      icp1b=ic-1
c
c  *** extract first word ***
c
      part1=runlb2(icp1a:icp1b)
      ilp1=str_length(part1)
c
c------------
c
      do ic=icp1b+1,ilnb
         if(runlb2(ic:ic).ne.' ') go to 30
      enddo
c
      return                            ! rest is blank
c
 30   continue
      icp2a=ic
c
      do ic=icp2a+1,ilnb
         if(runlb2(ic:ic).eq.' ') go to 40
      enddo
      ic=ilnb+1
c
 40   continue
      icp2b=ic-1
c
c  *** extract second word ***
c
      part2=runlb2(icp2a:icp2b)
      ilp2=str_length(part2)
c
c------------
c
      do ic=icp2b+1,ilnb
         if(runlb2(ic:ic).ne.' ') go to 50
      enddo
c
      return                            ! rest is blank
c
 50   continue
      icp3a=ic
c
C  *** last word(s) ***
c
      part3=runlb2(icp3a:ilnb)
      ilp3=str_length(part3)
c
      return
      end
 
