!
! This file collects all the routines that are used to perform a
! matrix-vector multiplication for COSMO and PCM. This includes
!
!   lx      : COSMO matrix
!   lstarx  : COSMO adjoint matrix
!
!   + service routines
!
!-------------------------------------------------------------------------------
! given a vector x, compute y = Lx, where L is the ddCOSMO matrix
! (off-diagonal blocks only).
!-------------------------------------------------------------------------------
!
subroutine lx(ddCosmo, n, x, y )
   !
   use ddcosmo_core , only : TddCosmo, zero, calcv, intrhs, prtsph
   !
   implicit none
   type(TddCosmo), intent(in) :: ddCosmo
   integer,                         intent(in)    :: n
   real*8,  dimension(ddCosmo%nylm, ddCosmo%nsph), intent(in)    :: x
   real*8,  dimension(ddCosmo%nylm, ddCosmo%nsph), intent(inout) :: y
   !
   integer             :: isph, istatus
   real*8, allocatable :: pot(:), vplm(:), basloc(:), vcos(:), vsin(:)
   !
   !-------------------------------------------------------------------------------
   !
   !     allocate workspaces
   allocate( pot(ddCosmo%ngrid), vplm(ddCosmo%nylm), basloc(ddCosmo%nylm), vcos(ddCosmo%lmax+1), &
      vsin(ddCosmo%lmax+1) , stat=istatus )
   if ( istatus.ne.0 ) then
      write(*, *) 'lx: allocation failed !'
      stop
   endif
   !
   if (ddCosmo%iprint.ge.5) call prtsph(ddCosmo, 'X', ddCosmo%nsph, 0, x)
   !
   !     initialize
   y = zero
   !
   !$omp parallel do default(shared) private(isph, pot, basloc, vplm, vcos, vsin) &
   !$omp schedule(dynamic)
   !
   !
   !     loop over spheres
   do isph = 1, ddCosmo%nsph
      !
      !       compute NEGATIVE action of off-digonal blocks
      call calcv(ddCosmo, .false., isph, pot, x, basloc, vplm, vcos, vsin)
      call intrhs(ddCosmo, isph, pot, y(:, isph))
      !
      !       action of off-diagonal blocks
      y(:, isph) = - y(:, isph)
      !
   enddo
   !
   if (ddCosmo%iprint.ge.5) call prtsph(ddCosmo, 'LX (off diagonal)', ddCosmo%nsph, 0, y)
   !
   !     deallocate workspaces
   deallocate( pot, basloc, vplm, vcos, vsin , stat=istatus )
   if ( istatus.ne.0 ) then
      write(*, *) 'lx: allocation failed !'
      stop
   endif
   !
   !
end subroutine lx
!-------------------------------------------------------------------------------
!
!
!
!
!
!-------------------------------------------------------------------------------
! given a vector x, compute y = L*x, where L* is the adjoint ddCOSMO matrix.
! if dodiag is set to .true., L includes the diagonal blocks, otherwise
! L only includes the off-diagonal ones.
!-------------------------------------------------------------------------------
!
subroutine lstarx(ddCosmo, n, x, y )
   !
   use ddcosmo_core , only : TddCosmo, zero, &
      adjrhs, prtsph
   !
   implicit none
   type(TddCosmo), intent(in) :: ddCosmo
   integer,                       intent(in)    :: n
   real*8,  dimension(ddCosmo%nylm, ddCosmo%nsph), intent(in)    :: x
   real*8,  dimension(ddCosmo%nylm, ddCosmo%nsph), intent(inout) :: y
   !
   integer             :: isph, ig, istatus
   real*8, allocatable :: xi(:, :), vplm(:), basloc(:), vcos(:), vsin(:)
   !
   !-------------------------------------------------------------------------------
   !
   !     allocate workspaces
   allocate( xi(ddCosmo%ngrid, ddCosmo%nsph), vplm(ddCosmo%nylm), basloc(ddCosmo%nylm), vcos(ddCosmo%lmax+1), &
      vsin(ddCosmo%lmax+1) , stat=istatus )
   if ( istatus.ne.0 ) then
      write(*, *) 'lstarx: allocation failed!'
      stop
   endif
   !
   if (ddCosmo%iprint.ge.5) call prtsph(ddCosmo, 'X', ddCosmo%nsph, 0, x)
   !
   !     initilize
   y = zero
   !
   !$omp parallel do default(shared) private(isph, ig)
   !
   !
   !     expand x over spherical harmonics
   !     ---------------------------------
   !
   !     loop over spheres
   do isph = 1, ddCosmo%nsph
      !
      !       loop over gridpoints
      do ig = 1, ddCosmo%ngrid
         !
         xi(ig, isph) = dot_product( x(:, isph), ddCosmo%basis(:, ig) )
         !
      enddo
   enddo
   !
   !$omp parallel do default(shared) private(isph, basloc, vplm, vcos, vsin) &
   !$omp schedule(dynamic)
   !
   !     compute action
   !     --------------
   !
   !     loop over spheres
   do isph = 1, ddCosmo%nsph
      !
      !       compute NEGATIVE action of off-digonal blocks
      call adjrhs(ddCosmo, isph, xi, y(:, isph), basloc, vplm, vcos, vsin)
      !
      !       action of off-diagonal blocks
      y(:, isph) = - y(:, isph)
      !
      !       add action of diagonal block
      !
   enddo
   !
   if (ddCosmo%iprint.ge.5) call prtsph(ddCosmo, 'L*X (off-diagonal)', ddCosmo%nsph, 0, y)
   !
   !     deallocate workspaces
   deallocate( xi, basloc, vplm, vcos, vsin , stat=istatus )
   if ( istatus.ne.0 ) then
      write(*, *) 'lstarx: allocation failed !'
      stop
   endif
   !
   !
end subroutine lstarx

!> given a vector x, apply the inverse diagonal (block) of the L matrix:
subroutine ldm1x(ddCosmo, n, x, y )
   use ddcosmo_core , only : TddCosmo
   implicit none
   type(TddCosmo), intent(in) :: ddCosmo
   integer, intent(in) :: n
   real*8, intent(in) :: x(ddCosmo%nylm, ddCosmo%nsph)
   real*8, intent(inout) :: y(ddCosmo%nylm, ddCosmo%nsph)
   integer                                        :: isph
   ! loop over spheres
   do isph = 1, ddCosmo%nsph
      ! apply inverse
      y(:, isph) = ddCosmo%facl*x(:, isph)
   enddo
end subroutine ldm1x
!-------------------------------------------------------------------------------
!
!
!
!
!
!-------------------------------------------------------------------------------
! compute the h^-1/2 norm of the increment on each sphere, then take the
! rms value.
!-------------------------------------------------------------------------------
!
real*8 function hnorm(ddCosmo, n, x )
   !
   use ddcosmo_core , only : TddCosmo, hsnorm
   !
   implicit none
   type(TddCosmo), intent(in) :: ddCosmo
   integer,                         intent(in) :: n
   real*8,  dimension(ddCosmo%nylm, ddCosmo%nsph), intent(in) :: x
   !
   integer                                     :: isph, istatus
   real*8                                      :: vrms, vmax
   real*8, allocatable                         :: u(:)
   !
   !-------------------------------------------------------------------------------
   !
   !     allocate workspace
   allocate( u(ddCosmo%nsph) , stat=istatus )
   if ( istatus.ne.0 ) then
      write(*, *) 'hnorm: allocation failed !'
      stop
   endif
   !
   !     loop over spheres
   do isph = 1, ddCosmo%nsph
      !
      !       compute norm contribution
      call hsnorm(ddCosmo, x(:, isph), u(isph) )
   enddo
   !
   !     compute rms of norms
   call rmsvec( ddCosmo%nsph, u, vrms, vmax )
   !
   !     return value
   hnorm = vrms
   !
   !     deallocate workspace
   deallocate( u , stat=istatus )
   if ( istatus.ne.0 ) then
      write(*, *) 'hnorm: deallocation failed !'
      stop
   endif
   !
   !
end function hnorm
!-------------------------------------------------------------------------------
