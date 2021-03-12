subroutine eqm_ball(ivec,zR,zZ,zphi,init,BR,BZ,BPHI,zpsi,ipsi,ierr)

  !  compute B vector on (R,Z) grid -- psi(R,Z) presumed known.

  !  **valid only for axisymmetric equilibria**

  !  axisymmetry is assumed... phi argument is ignored, but left in so
  !  that the routine can be used as an argument of the eqm_brz subroutine

  use xplasma_definitions
  use eqi_rzbox_module
  use i2mex_mod, ONLY : i2mex_o
  implicit NONE

  integer ivec                      ! vector lengths
  real*8 zR(ivec),zZ(ivec),zphi(ivec) ! locations to evaluate B vector

  integer init                      ! init-cleanup-compute control flag

  !  init=1  -- initialize
  !  init=2  -- cleanup
  !  init=0  -- compute

  !  ...a vector of B vectors...
  real*8 BR(ivec),BZ(ivec),BPHI(ivec) ! field vectors components (returned)

  real*8 zpsi(ivec)                 ! psi(poloidal) returned.
  !  for this routine, zpsi is a placeholder, it is not touched.
  !  the place holder is necessary to make this routine usable as a passed
  !  argument "userbvec" in xplasma_brz.

  integer ipsi                      ! =0 returned, no psi(poloidal) flag
  !  Psi is input, not output, in this routine.

  integer ierr                      ! completion code, 0=OK

  !--------------------------------

  real*8, dimension(:), allocatable :: zg,zrho,zchi,zRloc,zZloc,zdeti
  real*8, dimension(:,:), allocatable :: f2,f5
  real*8 zRaxis,zZaxis,zdum

  logical ioutsid(ivec),axisymm

  integer i,id_bmod,id_psi,id_g,id_Rgrid,id_Zgrid,id_R,id_Z,id_map

  integer :: id5(5),idth(5),idrho(5)
  integer :: id2(2),idR(2),idZ(2)
  integer,parameter :: idRdth=1,idRdrho=2,idZdth=3,idZdrho=4,idpsidrho=5

  integer :: nsnccwb,nsnccwi

  integer ivecin,ivecout
  integer, dimension(:), allocatable :: iregion

  real*8, parameter :: zrhomin = 1.0d-12
  !--------------------------------

  call xplasma_global_info(sp,ierr, axisymm=axisymm, &
       bphi_ccw=nsnccwb, jphi_ccw=nsnccwi)
  if(ierr.ne.0) return

  if(.not.axisymm) then
     ierr=107
     call xplasma_errmsg_append(sp,' ?eqm_ball: plasma must be axisymmetric!')
     go to 900
  endif

  ipsi=0

  !  cleanup call

  if(init.eq.2) then
     !  nothing need be done
     go to 900
  endif

  !  initialization; error checks done once...

  call xplasma_common_ids(sp,ierr, id_bmod=id_bmod,id_psi=id_psi, &
       id_g=id_g,id_R=id_R,id_Z=id_Z)

  if(init.eq.1) then

     if(min(id_bmod,id_psi,abs(nsnccwi),abs(nsnccwb)).eq.0) then
        call xplasma_errmsg_append(sp, &
             ' ?eqm_ball: core plasma equilibrium & fields must be defined!')
        ierr=608
        go to 900
     endif

     call xplasma_find_item(sp,'__RGRID',id_Rgrid,ierr,nf_noerr=.TRUE.)
     call xplasma_find_item(sp,'__ZGRID',id_Zgrid,ierr,nf_noerr=.TRUE.)

     if(min(id_Rgrid,id_Zgrid).eq.0) then
        call xplasma_errmsg_append(sp, &
             ' ?eqm_ball: R & Z external grids must be defined first.')
        ierr=615
        go to 900
     endif

     call xplasma_mag_axis(sp,ierr, zRaxis,zZaxis)
     if(ierr.ne.0) then
        ierr=9999
        call xplasma_errmsg_append(sp, &
             ' ?eqm_ball: unexpected error in xplasma_mag_axis call.')
        go to 900
     endif

     !  psi(R,Z) must exist

     call xplasma_eval_prof(sp,i2mex_o%id_psirz, &
          xplasma_R_coord,zRaxis,xplasma_Z_coord,zZaxis, &
          zdum,ierr)
     if(ierr.ne.0) then
        call xplasma_errmsg_append(sp, &
             ' ?eqm_ball: test evaluation of Psi(R,Z) failed.')
        go to 900
     endif

     !----------------------------------
     !  OK...
     go to 900                      ! exit initialization
  endif

  !  evaluation...

  allocate(zrho(ivec),zchi(ivec),zg(ivec),iregion(ivec))
  call eqi_fastinv(ivec,zR,zZ,zphi,zrhomin,zrho,zchi,iregion,ierr)
  if(ierr.ne.0) then
     deallocate(zrho,zchi,zg,iregion)
     go to 900
  endif

  do i=1,ivec
     if(zrho(i).gt.(1.0d0)) then
        !  outside core region
        ioutsid(i)=.TRUE.
        zrho(i)=1.00d0
     else
        !  inside core region and...
        ioutsid(i)=.FALSE.       ! core BR,BZ are available
     endif
  enddo

  !  ioutsid(i)=.FALSE. means (BR,BZ) will be evaluated by interpolating
  !    core BR(chi,phi) and BZ(chi,phi) functions;
  !  ioutsid(i)=.TRUE. means (BR,BZ) will be evaluated via 1/R * grad(psi)

  !  do Bphi now...

  call xplasma_eval_prof(sp,id_g,zrho,zg,ierr)
  if(ierr.ne.0) then
     deallocate(zrho,zchi,zg,iregion)
     go to 900
  endif

  do i=1,ivec
     Bphi(i)=nsnccwb*zg(i)/zR(i)
  enddo

  !  OK for (1/R)grad(Psi) ^ e_phi need to sort inside and outside pts...

  allocate(zRloc(ivec),zZloc(ivec))

  ivecin=0
  ivecout=0
  do i=1,ivec
     if(ioutsid(i)) then
        ivecout=ivecout+1
        zRloc(ivecout)=zR(i)
        zZloc(ivecout)=zZ(i)
     else
        ivecin=ivecin+1
        zrho(ivecin)=zrho(i)
        zchi(ivecin)=zchi(i)
     endif
  enddo
  PRINT *,ivecin,' pts inside; ',ivecout,' pts outside.'

  !  evaluate inside pts

  if(ivecin.gt.0) then

     !  set up to evaluate {R,Z,Psi} derivatives; settings must be consistent
     !  with parameters idRdrho, etc.
     id5(1:2)=id_R
     id5(3:4)=id_Z
     id5(5)=id_psi
     idth = (/1, 0, 1, 0, 0/)  ! d/dth controls
     idrho= (/0, 1, 0, 1, 1/)  ! d/drho controls

     allocate(f5(ivecin,5),zdeti(ivecin))

     call xplasma_eval_prof(sp,id5, &
          xplasma_theta_coord,zchi(1:ivecin),xplasma_rho_coord,zrho(1:ivecin),&
          f5,ierr, ideriv1s=idth, ideriv2s=idrho)

     if(ierr.eq.0) then
        ! fold in sign factor (direction of current) & 1/det[metric Jacobian]

        zdeti = nsnccwi/(f5(1:ivecin,idRdth)*f5(1:ivecin,idZdrho) - &
             f5(1:ivecin,idZdth)*f5(1:ivecin,idRdrho))

        !  copy back

        ivecin=0
        do i=1,ivec
           if(.not.ioutsid(i)) then
              ivecin=ivecin+1
              BR(i)=f5(ivecin,idpsidrho)*f5(ivecin,idRdth)*zdeti(ivecin)/zR(i)
              BZ(i)=f5(ivecin,idpsidrho)*f5(ivecin,idZdth)*zdeti(ivecin)/zR(i)
           endif
        enddo
     endif
     deallocate(f5,zdeti)
  endif
  if(ierr.ne.0) then
     deallocate(zrho,zchi,zg,iregion,zRloc,zZloc)
     go to 900
  endif

  !  evaluate outside points using gradpsi & note sign factor...

  if(ivecout.gt.0) then

     !  lookup on psi's R,Z grid

     allocate(f2(ivecout,2))

     id2(1:2)=i2mex_o%id_psirz
     idR = (/ 1, 0 /)
     idZ = (/ 0, 1 /)

     call xplasma_eval_prof(sp,id2, &
          xplasma_R_coord,zRloc(1:ivecout),xplasma_Z_coord,zZloc(1:ivecout), &
          f2,ierr, ideriv1s=idR, ideriv2s=idZ)

     if(ierr.eq.0) then

        !  f2(:,1) = dPsi/dR
        !  f2(:,2) = dPsi/dZ

        ivecout=0
        do i=1,ivec
           if(ioutsid(i)) then
              ivecout=ivecout+1
              BR(i)=nsnccwi*f2(ivecout,2)/zRloc(ivecout)  ! or /zR(i)
              BZ(i)=-nsnccwi*f2(ivecout,1)/zRloc(ivecout) ! or /zR(i)
           endif
        enddo
     endif

  endif

  deallocate(zrho,zchi,zg,iregion,zRloc,zZloc)
  if(ierr.ne.0) then
     go to 900
  else
     go to 1000
  endif

900 continue
  bphi=0.0d0
  bz=0.0d0
  br=0.0d0

1000 continue

end subroutine eqm_ball
