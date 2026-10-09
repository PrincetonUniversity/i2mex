subroutine phcdf_define_int(ncid,name,inval)
  !  phcdf_define -- Define Variables and Dimensions for NetCDF Dataset
  !    routine is called from generated subroutine trdatbuf_def_cdf
  !   Routine supports only scalar and 1D array integer data
  !
  implicit none
  include 'netcdf.inc'
 
  integer ncid, inval, nvdims, dimid, status, varid
  character*(*) name
  integer, dimension(1) :: vdims ! vector of dimension IDs corresponding to variable dimensions
  character*13 :: dimnam
  integer :: lunmsg_tdb ! integer function lunmsg_tdb(idum)
  !
  !  ncid -- identifier of open ph.cdf file
  !  name -- character*(*) name of item to write
  !  inval -- dimension  of item (integer array) to write
 
  if(inval.eq.1)then ! define integer scalar
    nvdims = 0
  elseif(inval.gt.1)then ! define 1D integer array
    nvdims = 1
    !
    !  Create and Define Dimension
    if(inval.le.99999)then
      write(dimnam,"('dim_',i5.5)") inval
    elseif(inval.le.999999999)then
      write(dimnam,"('dim_',i9.9)") inval
    else
      write(lunmsg_tdb(0),*) ' ?phcdf_define_int: dimension > 1e9 not supported '
      return
    endif
    !
    status = nf_inq_dimid(ncid, dimnam, dimid)
    if(status .ne. NF_NOERR)then   ! define dimension
       status = nf_def_dim(ncid, dimnam, inval, dimid)
    endif
    vdims(1) = dimid  ! dimension ID for first (and only) dimension
  else
      write(lunmsg_tdb(0),*) ' ?phcdf_define_int: dimension < 0 not supported '
      return
  endif
  status = nf_def_var(ncid, name, NF_INT, nvdims, vdims, varid)
  !
  return
end subroutine phcdf_define_int

subroutine phcdf_write_int(ncid,name,ival,inval)
  !  phcdf_write -- write operation on ph.cdf file
  !    routine is called from generated subroutine trdatbuf_wr_cdf
  !
  implicit none
  include 'netcdf.inc'
 
  integer inval,ncid, status, varid
  character*(*) name
  integer ival(inval)
  !
  !  ncid -- identifier of open ph.cdf file
  !  name -- character*(*) name of item to write
  !  ival(inval) -- dimension and value of item (integer array) to write
 
  status = nf_inq_varid(ncid, name, varid)
  status = nf_put_var_int(ncid, varid, ival)
  !
  return
end subroutine phcdf_write_int

subroutine phcdf_write_int1(ncid,name,ival)
  !  phcdf1_write -- write operation on ph.cdf file
  !    routine is called from generated subroutine trdatbuf_wr_cdf
  !
  implicit none

  integer, intent(in) :: ncid
  character*(*), intent(in) :: name
  integer, intent(in) :: ival
  integer ival2(1)
  ival2(1) = ival
  call phcdf_write_int(ncid,name,ival2,1)
  !
  return
end subroutine phcdf_write_int1

subroutine phcdf_read_int(ncid,name,ival,inval,ier)
  !  phcdf_read -- read operation on ph.cdf file
  !    routine is called from generated subroutine trdatbuf_rd_cdf
  !
  implicit none
  include 'netcdf.inc'
 
  integer inval,ier,ncid, status, varid
  character*(*) name
  integer ival(inval)
  !
  !  ncid -- identifier of open ph.cdf file
  !  name -- character*(*) name of item to read
  !  ival(inval) -- dimension and value of item (integer array) to read
  !  ier -- error code 0=OK, 1=Data Not Found, 2=Error Reading Data
 
  ier = 0
  ival = 0 ! clear item to zero
  status = nf_inq_varid(ncid, name, varid)
  if(status.eq.NF_NOERR)then ! data present, now read it
    status = nf_get_var_int(ncid, varid, ival)
    if(status.ne.NF_NOERR)ier = 2 ! error reading data
  else
    ier = 1 ! no data found
    return
  endif
  !
  return
end subroutine phcdf_read_int

subroutine phcdf_read_int1(ncid,name,ival,ier)
  !  phcdf1_read -- read operation on ph.cdf file
  !    routine is called from generated subroutine trdatbuf_rd_cdf
  !
  implicit none

  integer, intent(in) :: ier, ncid
  character*(*), intent(in) :: name
  integer, intent(out) :: ival
  integer ival2(1)
  !
  call phcdf_read_int(ncid,name,ival2,1,ier)
  ival = ival2(1)
  !
  return
end subroutine phcdf_read_int1

subroutine phcdf_define_double(ncid,name,inval)
  !  phcdf_define -- Define Variables and Dimensions for NetCDF Dataset
  !    routine is called from generated subroutine trdatbuf_def_cdf
  !   Routine supports only scalar and 1D array double data
  !
  implicit none
  include 'netcdf.inc'
 
  integer ncid, inval, nvdims, dimid, status, varid
  character*(*) name
  integer, dimension(1) :: vdims ! vector of dimension IDs corresponding to variable dimensions
  character*13 :: dimnam
     integer :: lunmsg_tdb ! integer function lunmsg_tdb(idum)
  !
  !  ncid -- identifier of open ph.cdf file
  !  name -- character*(*) name of item to write
  !  inval -- dimension  of item (real*8 array) to write
 
  if(inval.eq.1)then ! define real*8 scalar
    nvdims = 0
  elseif(inval.gt.1)then ! define 1D real*8 array
    nvdims = 1
    !
    !  Create and Define Dimension
    if(inval.le.99999)then
      write(dimnam,"('dim_',i5.5)") inval
    elseif(inval.le.999999999)then
      write(dimnam,"('dim_',i9.9)") inval
    else
      write(lunmsg_tdb(0),*) ' ?phcdf_define_double: dimension > 1e9 not supported '
      return
    endif
    !
    status = nf_inq_dimid(ncid, dimnam, dimid)
    if(status .ne. NF_NOERR)then   ! define dimension
       status = nf_def_dim(ncid, dimnam, inval, dimid)
    endif
    vdims(1) = dimid  ! dimension ID for first (and only) dimension
  else
      write(lunmsg_tdb(0),*) ' ?phcdf_define_double: dimension < 0 not supported '
      return
  endif
  status = nf_def_var(ncid, name, NF_DOUBLE, nvdims, vdims, varid)
  !
  return
end subroutine phcdf_define_double
 
subroutine phcdf_write_double(ncid,name,val,inval)
  !  phcdf_write -- write operation on ph.cdf file
  !    routine is called from generated subroutine trdatbuf_wr_cdf
  !
  implicit none
  include 'netcdf.inc'
 
  integer inval,ncid, status, varid
  character*(*) name
  real*8 val(inval)
  !
  !  ncid -- identifier of open ph.cdf file
  !  name -- character*(*) name of item to write
  !  val(inval) -- dimension and value of item (real*8 array) to write
 
  status = nf_inq_varid(ncid, name, varid)
  status = nf_put_var_double(ncid, varid, val)
  !
  return
end subroutine phcdf_write_double
 
subroutine phcdf_read_double(ncid,name,val,inval,ier)
  !  phcdf_read -- read operation on ph.cdf file
  !    routine is called from generated subroutine trdatbuf_rd_cdf
  !
  implicit none
  include 'netcdf.inc'
 
  integer inval,ier,ncid, status, varid
  character*(*) name
  real*8 val(inval)
  !
  !  ncid -- identifier of open ph.cdf file
  !  name -- character*(*) name of item to read
  !  val(inval) -- dimension and value of item (real*8 array) to read
  !  ier -- error code 0=OK, 1=Data Not Found, 2=Error Reading Data
 
  ier = 0
  status = nf_inq_varid(ncid, name, varid)
  if(status.eq.NF_NOERR)then ! data present, now read it
    status = nf_get_var_double(ncid, varid, val)
    if(status.ne.NF_NOERR)ier = 2 ! error reading data
  else
    ier = 1 ! no data found
    return
  endif
  !
  return
end subroutine phcdf_read_double
subroutine trdatbuf_open(filename,mode,ncid)
  !
  ! OPEN trdatbuf NetCDF file
  !
  implicit none
  include 'netcdf.inc'
  character*(*), intent(in) :: filename ! path to NetCDF file
  character*(*), intent(in) :: mode     ! "r" or "R" for read
                                        ! "w" or "W" for write
  integer, intent(out) :: ncid          ! NetCDF channel # returned
  ! note: if the NetCDF open fails, ncid=0 is returned.

  integer :: status
 
  ncid=0
  if(mode.eq.'w' .or.  mode.eq.'W')then
    status = nf_create(filename, NF_CLOBBER, ncid)
    if(status.ne.NF_NOERR) ncid=0
  elseif(mode.eq.'r' .or.  mode.eq.'R')then
    status = nf_open(filename, NF_NOWRITE, ncid)
    if(status.ne.NF_NOERR) ncid=0
  endif
 
end subroutine trdatbuf_open
 
!-----------------------------------------
 
subroutine trdatbuf_close(ncid)
  !
  ! CLOSE trdatbuf NetCDF file
  !
  implicit none
  include 'netcdf.inc'
  integer, intent(inout) :: ncid

  integer :: status
 
  status = nf_close(ncid)
  ncid=0
 
end subroutine trdatbuf_close
