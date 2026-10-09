C-----------------------------------------------------------------------
C  gettrng -- find range of pts btw 2 values
C
      subroutine gettrng(ztime1,ztime2,tarray,narray,it1,it2)
C
      real ztime1,ztime2,tarray(narray)
      integer it1,it2
C
C  assume ztime2 > ztime1
C  assume tarray already sorted in increasing order
C
C---------------------------
C
      it1=0
      it2=0
      do it=1,narray
         if(tarray(it).ge.ztime1) go to 10
      enddo
C
      return       ! nothing found.
C
 10   continue
      if(tarray(it).gt.ztime2) return    ! nothing found.
C
      it1=it
      do it=it1,narray
         if(tarray(it).gt.ztime2) go to 20
      enddo
      it=narray+1
C
 20   continue
      it2=it-1
C
      return
      end
