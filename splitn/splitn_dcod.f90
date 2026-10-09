subroutine splitn_rdcod(word,rbuf,num,ier)

  ! REAL decode of character string

  implicit NONE

  character*(*), intent(in) :: word   ! string to decode: no leading or
  !                                   ! trailing blanks

  integer, intent(in) :: num          ! no. of copies to store (repeat count)
  real, intent(out) :: rbuf(num)      ! where to store results

  integer, intent(out) :: ier         ! completion code (0=OK)

  !------------------------------------------
  ! if an error occurs, rbuf is not modified and a message
  ! is written on unit 6.
  ! a DECIMAL POINT must be contained in the input string.
  !------------------------------------------

  integer ildp,ilen

  character*30 jbuf

  real zval

  !------------------------------------------

  ier=0
  ildp = index(word,'.')
  if(ildp.le.0) then
     write(6,*) ' ?splitn_dcod: no decimal point in floating point number: "',&
          word,'"'
     ier=2
     return
  endif

  ilen=len(word)
  if(ilen.gt.30) then
     write(6,*) ' ?splitn_dcod: floating point number too long: ',word
     ier=3
     return
  endif

  jbuf=' '
  jbuf(30-ilen+1:30)=word
  call uupper(jbuf)
  if(jbuf(28:30).eq.'_R8') then
     read(jbuf(1:27),'(G27.0)',iostat=ier) zval
  else
     read(jbuf,'(G30.0)',iostat=ier) zval
  endif

  if(ier.ne.0) then
     write(6,*) ' ?splitn_dcod: error decoding as floating point: "', &
          word,'"'
     ier=1
     return
  endif

  ! OK

  rbuf(1:num)=zval

end subroutine splitn_rdcod

subroutine splitn_ddcod(word,dbuf,num,ier)

  ! REAL*8 decode of character string

  implicit NONE

  character*(*), intent(in) :: word   ! string to decode: no leading or
  !                                   ! trailing blanks

  integer, intent(in) :: num          ! no. of copies to store (repeat count)
  real*8, intent(out) :: dbuf(num)    ! where to store results

  integer, intent(out) :: ier         ! completion code (0=OK)

  !------------------------------------------
  ! if an error occurs, dbuf is not modified and a message
  ! is written on unit 6.
  ! a DECIMAL POINT must be contained in the input string.
  !------------------------------------------

  integer ildp,ilen

  character*30 jbuf

  real*8 zval

  !------------------------------------------

  ier=0
  ildp = index(word,'.')
  if(ildp.le.0) then
     write(6,*) ' ?splitn_dcod: no decimal point in floating point number: "',&
          word,'"'
     ier=2
     return
  endif

  ilen=len(word)
  if(ilen.gt.30) then
     write(6,*) ' ?splitn_dcod: floating point number too long: ',word
     ier=3
     return
  endif

  jbuf=' '
  jbuf(30-ilen+1:30)=word
  call uupper(jbuf)
  if(jbuf(28:30).eq.'_R8') then
     read(jbuf(1:27),'(G27.0)',iostat=ier) zval
  else
     read(jbuf,'(G30.0)',iostat=ier) zval
  endif

  if(ier.ne.0) then
     write(6,*) ' ?splitn_dcod: error decoding as floating point: "', &
          word,'"'
     ier=1
     return
  endif

  ! OK

  dbuf(1:num)=zval

end subroutine splitn_ddcod

subroutine splitn_idcod(word,ibuf,num,ier)

  ! INTEGER decode of character string

  implicit NONE

  character*(*), intent(in) :: word   ! string to decode: no leading or
  !                                   ! trailing blanks

  integer, intent(in) :: num          ! no. of copies to store (repeat count)
  integer, intent(out) :: ibuf(num)   ! where to store results

  integer, intent(out) :: ier         ! completion code (0=OK)

  !------------------------------------------
  ! if an error occurs, ibuf is not modified and a message
  ! is written on unit 6.
  !------------------------------------------

  integer ilen

  character*30 jbuf

  integer ival

  !------------------------------------------

  ier=0
  ilen=len(word)
  if(ilen.gt.30) then
     write(6,*) ' ?splitn_dcod: integer too long: ',word
     ier=3
     return
  endif

  jbuf=' '
  jbuf(30-ilen+1:30)=word
  read(jbuf,'(I30)',iostat=ier) ival

  if(ier.ne.0) then
     write(6,*) ' ?splitn_dcod: error decoding as integer: "', &
          word,'"'
     ier=1
     return
  endif

  ! OK

  ibuf(1:num)=ival

end subroutine splitn_idcod

subroutine splitn_ldcod(word,ibuf,num,ier)

  ! LOGICAL decode of character string

  implicit NONE

  character*(*), intent(in) :: word   ! string to decode: no leading or
  !                                   ! trailing blanks

  integer, intent(in) :: num          ! no. of copies to store (repeat count)
  logical, intent(out) :: ibuf(num)   ! where to store results

  integer, intent(out) :: ier         ! completion code (0=OK)

  !------------------------------------------
  ! if an error occurs, ibuf is not modified and a message
  ! is written on unit 6.
  !------------------------------------------

  integer ilen

  character*30 jbuf

  logical ival

  !------------------------------------------

  ier=0
  ilen=len(word)
  if(ilen.gt.30) then
     write(6,*) ' ?splitn_dcod: logical value too long: ',word
     ier=3
     return
  endif

  jbuf=' '
  jbuf(30-ilen+1:30)=word
  read(jbuf,'(L30)',iostat=ier) ival

  if(ier.ne.0) then
     write(6,*) ' ?splitn_dcod: error decoding as logical(T/F): "', &
          word,'"'
     ier=1
     return
  endif

  ! OK

  ibuf(1:num)=ival

end subroutine splitn_ldcod

subroutine splitn_cdcod(word,maxlen,cbuf,num,ier)

  ! CHARACTER decode of character string

  implicit NONE

  character*(*), intent(in) :: word   ! string to decode: no leading or
  !                                   ! trailing blanks; enclosed in quotes!

  integer, intent(in) :: maxlen       ! maximum length of enclosed string
  !                                   ! not counting trailing blanks

  integer, intent(in) :: num          ! no. of copies to store (repeat count)
  character*(*), intent(out) :: cbuf(num)   ! where to store results

  integer, intent(out) :: ier         ! completion code (0=OK)

  !------------------------------------------
  ! if an error occurs, cbuf is not modified and a message
  ! is written on unit 6.
  !------------------------------------------

  integer ilen
  character*128 cval

  !------------------------------------------

  ier=0
  ilen=len(word)

  cval=' '
  if( ((word(1:1).ne."'").and.(word(1:1).ne.'"')) .or. &
       (word(1:1).ne.word(ilen:ilen)) ) then
     write(6,*) ' ?splitn_dcod: string not enclosed in matching quotes: ', &
          word
     ier=1
     return
  endif

  ! OK

  if(ilen.le.len(cval)+2) then
     cval=word(2:ilen-1)
     ilen=ilen-2
  else
     write(6,*) " ?splitn_dcod: warning: internal buffer too small."
     ilen=maxlen+1
  endif

  if(ilen.gt.maxlen) then
     write(6,*) " ?splitn_dcod: string length error."
     ier=2
  else
     ier=0
  endif

  cbuf(1:num)=cval

end subroutine splitn_cdcod
