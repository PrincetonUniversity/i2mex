      subroutine fwc_delete1(path_expr,ier)
C
C  access to fwc_delete, with single argument of form <path>/<wc-expr>
C  break the string argument into two pieces; "*" must be in the filename
C  tail; no support yet for wildcard in directory path
C
      implicit NONE
C
      character*(*), intent(in) :: path_expr
C
      integer, intent(out) :: ier
C
C---------------------
      integer :: ilen,ic,ipath_end
C---------------------
      ier = 0

      ilen = len(trim(path_expr))
      if(ilen.eq.0) then
         write(6,*) ' ?fwc_delete1 -- path_expr argument is blank.'
         ier = 1
         return
      endif
C
      ipath_end=0
      do ic=ilen,1,-1
         if(path_expr(ic:ic).eq.'/') then
            ipath_end=ic
            exit
         endif
      enddo
C
      if(ipath_end.eq.ilen) then
         call fwc_delete(path_expr(1:ilen),' ',ier)
      else if(ipath_end.eq.0) then
         call fwc_delete(' ',path_expr(1:ilen),ier)
      else
         call fwc_delete(path_expr(1:ipath_end),
     >        path_expr(ipath_end+1:ilen),ier)
      endif
C
      return
      end
C------------------------------------------------------------------------
      subroutine fwc_delete(path,expr,ier)
C
C  dmc 16 June 2011 -- clean out selected files in a directory
C    (hidden files & subdirectories not included)
C
C    expr = * -- remove all files
C    expr = <prefix>* -- remove all files starting w/string <prefix>
C    expr = *<suffix> -- remove all files ending w/string <suffix>
C    expr = <prefix>*<suffix> -- remove all files starting & ending
C           w/strings <prefix> and <suffix> respectively
C
C  wild card symbol "*" must appear exactly once in expr
C
      use iso_c_binding, only: c_char, c_null_char
      implicit NONE
C
      character*(*), intent(in) :: path  ! directory to be cleaned...
      character*(*), intent(in) :: expr  ! wildcard expression
      integer, intent(out) :: ier        ! completion status, 0=normal
C
C  if (path) is blank, remove files in $cwd !!
C
C  return error code
C    ier=0 -- success
C    ier=1 -- error
C
      integer ilp,ilw,istar,nstar,ii
C
      character*200 :: zexpr
C
      character(kind=c_char) cpath(500),cprefix(100),csuffix(100)
      integer cclean_dir,cwc_delete
C
C----------------------------------
C
      ier=0
C
      ilp=len_trim(path)
C
      if(index(path(1:ilp),'*').gt.0) then
         write(6,*) ' ?fwc_delete:  wildcard in directory path: ',
     >        ' not supported.'
         write(6,*) '  '//path(1:ilp)
         ier = 1
         return
      endif
C
      if(len(adjustl(trim(expr))).gt.len(zexpr)) then
         write(6,*) ' ?fwc_delete:  wildcard expression too long: ',
     >        trim(expr)
         ier=1
         return
      endif
C
      zexpr = adjustl(trim(expr))
C
      ilw=len_trim(zexpr)
      if(ilw.eq.0) then
         write(6,*) ' ?fwc_delete:  wildcard expression is blank.'
         ier=1
         return
      endif
C
      nstar=0
      do ii=1,ilw
         if(zexpr(ii:ii).eq.'*') then
            nstar=nstar+1
            istar=ii
         endif
      enddo
C
      if(nstar.eq.0) then
         write(6,*) ' ?fwc_delete: "*" missing in wildcard expression:',
     >        ' '//trim(zexpr)
         ier=1
         return
      endif
C
      if(nstar.gt.1) then
         write(6,*) ' ?fwc_delete: multiple "*" symbols in ',
     >        'wildcard expression: '//trim(zexpr)
         ier=1
         return
      endif
C
      if(ilp.eq.0) then
         call cstring('./',cpath,'2C')
      else if(path(ilp:ilp).ne.'/') then
         call cstring(path(1:ilp)//'/',cpath,'2C')
      else
         call cstring(path(1:ilp),cpath,'2C')
      endif
C
      if(ilw.eq.1) then
         !  bare "*": use cclean_dir(cpath), delete all files...
         ier=cclean_dir(cpath)
      else
C
C  prefix is zexpr(1:istar-1); suffix is zexpr(istar+1:ilw)
C  (can be null)
C
         cprefix(1)=c_null_char
         csuffix(1)=c_null_char
         if(istar.gt.1) then
            call cstring(zexpr(1:istar-1),cprefix,'2C')
         endif
         if(istar.lt.ilw) then
            call cstring(zexpr(istar+1:ilw),csuffix,'2C')
         endif
C
         ier=cwc_delete(cpath,cprefix,csuffix)
      endif
C
      return
      end
