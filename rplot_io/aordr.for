C******************** START FILE AORDR.FOR ; GROUP PLOTR1 ******************
C-----------------------------------------------------------------
C  AORDR -- CALCULATE ALPHABETIC ORDER ON SERIES OF 10+ CHARACTER
C   ABREVIATIONS -- no duplication check
C
C   side effect -- ALL ABBREVIATIONS CAPITALIZED!
C
      SUBROUTINE AORDR(IORD,ABB,N)
      INTEGER IORD(N)
      CHARACTER*(*) ABB(N)
C
      IF(N.LE.0) RETURN
      do i=1,n
         call uupper(abb(i))
      enddo
C
      IORD(1)=1
C
      IF(N.EQ.1) RETURN
C
      DO 500 INEX=2,N
         call aordr_add(abb,iord,inex)
 500  CONTINUE
C
      RETURN
      END
C-------------------------
      subroutine aordr_add(abb,iord,inex)
C
C  given ordering of elements 1...(inex-1), add in inex'th element.
C
      character*(*) abb(inex)
      integer iord(inex)
C
      IN=INEX-1
      if(in.eq.0) then
         imid=1
         go to 500
      endif
C  BINARY SEARCH FOR PLACE IN ORDER
      IMID=I_AORDR(ABB,IORD,IN,ABB(INEX))
      if(imid.eq.inex) go to 500
C
      if(abb(iord(imid)).eq.abb(inex)) then
         write(lunzer(0),1999) abb(inex)
 1999    format(' %aordr_add:  warning:  duplicate name:  ',a)
      endif
C
      DO 300 IL=IN,IMID,-1
         ILP1=IL+1
         IORD(ILP1)=IORD(IL)
 300  CONTINUE
 500  CONTINUE
      IORD(IMID)=INEX
C
      RETURN
C
      END
C-------------------------
      subroutine aordr_del(abb,iord,indim,inum,abbdel)
C
C  remove indicated item from list & ordering; decrement list length
C
      character*(*) abb(indim)          ! list array
      integer iord(indim)               ! ordering array
      integer inum                      ! actual length of list (decremented)
      character*(*) abbdel              ! item to delete
C
      iloc=i_aordr(abb,iord,inum,abbdel)
      if(abb(iord(iloc)).ne.abbdel) then
         write(lunzer(0),1999) abbdel
 1999    format(' %aordr_del:  warning:  no such name:  ',a)
         return
      endif
C
      ipos=iord(iloc)
      do i=iloc,inum-1
         iord(i)=iord(i+1)
      enddo
      do i=ipos,inum-1
         abb(i)=abb(i+1)
      enddo
C
      inum=inum-1
      do i=1,inum
         if(iord(i).gt.ipos) iord(i)=iord(i)-1
      enddo
C
      return
      end
C-------------------------
      integer function i_aordr(abb,iord,in,abnew)
C
C  function returns position in list, into which abnew should be
C  inserted.  return in+1 if abnew is .gt. last item in abb
C  otherwise return address of element before which abnew should
C  be inserted; no duplication check; abnew equal to indicated
C  element is possible!
C
      character*(*) abb(in)             ! sorted list
      integer iord(in)                  ! sorted list ordering
      character*(*) abnew               ! new item for list
C
      ITOP=1
      IBOT=IN
      IMID=IN+1
      IF(ABNEW.GT.ABB(IORD(IN))) GO TO 200
 20   CONTINUE
      IMID=(ITOP+IBOT)/2
C
      IF(ABNEW.GT.ABB(IORD(IMID))) GO TO 100
C  ABB COMES HIGHER IN LIST
      IBOT=IMID
C  SEARCH COMPLETE?
      IF(IBOT.EQ.ITOP) GO TO 200
      GO TO 20                          ! NOPE
C  ABB COMES FURTHER DOWN IN LIST
 100  CONTINUE
      ITOP=IMID+1
      GO TO 20
C  ORDER-- LOCATION IMID
 200  CONTINUE
      i_aordr=IMID
      return
      end
C
      integer function ifind_ordr(abb,iord,inum,abbrev)
C
      character*(*) abb(inum)           ! the list
      integer iord(inum)                ! the ordering
      character*(*) abbrev              ! the function sought
C
C  find function id using sort array; look for exact match
C  return 0 if no match
C
      if(inum.eq.0) then
         ifind_ordr=0
      else
         ii=i_aordr(abb,iord,inum,abbrev)
         if((ii.le.0).or.(ii.gt.inum)) then
            ifind_ordr=0
         else if(abb(iord(ii)).ne.abbrev) then
            ifind_ordr=0
         else
            ifind_ordr=iord(ii)
         endif
      endif
C
      return
      end
C******************** END FILE AORDR.FOR ; GROUP PLOTR1 ******************
