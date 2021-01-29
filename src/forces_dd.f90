subroutine forces_dd(ddCosmo, n, phi, sigma, s, fx)
   use ddcosmo_core
   !
   ! Sample driver for the calculation of the ddCOSMO forces.                     !
   !
   implicit none
   !
   type(TddCosmo), intent(in) :: ddCosmo
   integer,                         intent(in)    :: n
   real*8,  dimension(ddCosmo%ncav),        intent(in)    :: phi
   real*8,  dimension(ddCosmo%nylm, ddCosmo%nsph),   intent(in)    :: sigma, s
   real*8,  dimension(3, n),         intent(inout) :: fx
   !
   integer :: isph, ig, ii, c1, c2, cr
   real*8  :: fep
   !
   real*8, allocatable :: xi(:, :), phiexp(:, :), zeta(:), ef(:, :)
   real*8, allocatable :: basloc(:), dbsloc(:, :), vplm(:), vcos(:), vsin(:)
   !
   allocate (xi(ddCosmo%ngrid, ddCosmo%nsph), phiexp(ddCosmo%ngrid, ddCosmo%nsph))
   allocate (basloc(ddCosmo%nylm), dbsloc(3, ddCosmo%nylm), vplm(ddCosmo%nylm), vcos(ddCosmo%lmax+1), vsin(ddCosmo%lmax+1))
   !
   ! initialize the timer:
   !
   call system_clock(count_rate=cr)
   call system_clock(count=c1)
   !
   ! compute xi:
   !
   !$omp parallel do default(shared) private(isph, ig)
   do isph = 1, ddCosmo%nsph
      do ig = 1, ddCosmo%ngrid
         xi(ig, isph) = dot_product(s(:, isph), ddCosmo%basis(:, ig))
      end do
   end do
   !$omp end parallel do
   !
   if (ddCosmo%iprint.ge.4) call ptcart(ddCosmo, 'xi', ddCosmo%nsph, 0, xi)
   !
   ! expand the potential on a sphere-by-sphere basis (needed for parallelism):
   !
   ii = 0
   phiexp = zero
   do isph = 1, ddCosmo%nsph
      do ig = 1, ddCosmo%ngrid
         if (ddCosmo%ui(ig, isph).gt.zero) then
            ii = ii + 1
            phiexp(ig, isph) = phi(ii)
         end if
      end do
   end do
   !
   fx = zero
   do isph = 1, ddCosmo%nsph
      call fdoka(ddCosmo, isph, sigma, xi(:, isph), basloc, dbsloc, vplm, vcos, vsin, fx(:, isph))
      call fdokb(ddCosmo, isph, sigma, xi, basloc, dbsloc, vplm, vcos, vsin, fx(:, isph))
      call fdoga(ddCosmo, isph, xi, phiexp, fx(:, isph))
   end do
   !
   2000 format(1x, 'ddCOSMO-only contributions to the forces (atomic units):', /, &
      1x, ' atom', 15x, 'x', 15x, 'y', 15x, 'z')
   !
   if (ddCosmo%iprint.ge.4) then
      write(iout, 2000)
      do isph = 1, ddCosmo%nsph
         write(6, '(1x, i5, 3f16.8)') isph, fx(:, isph)
      end do
   end if
   !
   deallocate (basloc, dbsloc, vplm, vcos, vsin)
   !
   call system_clock(count=c2)
   if (ddCosmo%iprint.gt.0) then
      write(iout, 1010) dble(c2-c1)/dble(cr)
      1010 format(' the computation of the ddCOSMO part of the forces took ', f8.3, ' seconds.')
   end if
   !
   deallocate (xi, phiexp)
   !
   ! scale the forces time the cosmo factor:
   !
   fep = pt5*(ddCosmo%eps-one)/ddCosmo%eps
   fx  = fep*fx
   !
   return
   end
