program main
   use ddcosmo_core
   !                                                                              !
   ! Sample interface with the ddcosmo module.                                    !
   ! This program reads the various ddcosmo parameters and the geometrical        !
   ! information of the molecule from a text file. The solute is described        !
   ! by point charges, which are also read from a file.                           !
   ! Finally, the van der Waals radii are also read in.                           !
   !                                                                              !
   ! This program shows how to call the ddcosmo routines and how to initialize    !
   ! the various ddcosmo quantities; also, shows where the user should place      !
   ! his/her routines to compute the molecular electrostatic quantities.          !
   !                                                                              !
   implicit none
   !
   character(len=:), allocatable :: arg
   type(TddCosmo) :: ddCosmo
   integer :: i, ii, isph, ig, n
   real*8  :: tobohr, esolv, xx(1)
   real*8, parameter :: toang=0.52917721092d0, tokcal=627.509469d0
   !
   ! quantities to be allocated by the user.
   ! - solute's parameters, such as coordinates, vdw radii and
   !   point charges used to model the solute (or multipoles, or
   !   qm density...)
   !
   real*8, allocatable :: x(:), y(:), z(:), rvdw(:), charge(:)
   !
   ! - electrostatic potential phi(ncav) and psi vector psi(nylm, n)
   !
   real*8, allocatable :: phi(:), psi(:, :)
   !
   ! - ddcosmo solution sigma (nylm, n) and adjoint solution s(nylm, n)
   !
   real*8, allocatable :: sigma(:, :), s(:, :)
   !
   ! - forces:
   !
   real*8, allocatable :: fx(:, :), zeta(:), ef(:, :)
   !
   ! - for qm solutes, fock matrix contribution.
   !
   if (command_argument_count() /= 1) then
      error stop "Please provide one input file as command line argument"
   end if
   call get_argument(1, arg)
   ! here, we read all the ddcosmo parameters from a file named Input.txt
   !
   open (unit=100, file=arg, form='formatted', access='sequential')
   !
   ! scalar parameters. the variables are defined in the ddcosmo module and are common to
   ! all the ddcosmo routines (no need to declare them if ddcosmo.mod is loaded.)
   !
   read(100, *) ddCosmo%iprint      ! printing flag
   read(100, *) ddCosmo%nproc       ! number of openmp threads
   read(100, *) ddCosmo%lmax        ! max angular momentum of spherical harmonics basis
   read(100, *) ddCosmo%ngrid       ! number of lebedev points
   read(100, *) ddCosmo%iconv       ! 10^(-iconv) is the convergence threshold for the iterative solver
   read(100, *) ddCosmo%igrad       ! whether to compute (1) or not (0) forces
   read(100, *) ddCosmo%eps         ! dielectric constant of the solvent
   read(100, *) ddCosmo%eta         ! regularization parameter
   !
   read(100, *) n           ! number of atoms
   !
   allocate (x(n), y(n), z(n), rvdw(n), charge(n))
   !
   ! we also read from the same file the charges, coordinates and vdw radii.
   ! in this example, the coordinates and radii are read in angstrom and
   ! converted in bohr before calling ddinit.
   !
   do i = 1, n
      read(100, *) charge(i), x(i), y(i), z(i), rvdw(i)
   end do
   tobohr = 1.0d0/toang
   x    = x*tobohr
   y    = y*tobohr
   z    = z*tobohr
   rvdw = rvdw*tobohr
   !
   close (100)
   !
   ! call the initialization routine. this routine allocates memory, computes some
   ! quantities for internal use and creates and fills an array ccav(3, ncav) with
   ! the coordinates of the grid points at which the user needs to compute the potential.
   ! ncav is the number of external grid points and nylm the number of spherical
   ! harmonics functions used for the expansion of the various ddcosmo quantities;
   ! both are computed by ddinit and defined as common variables in ddcosmo.mod.
   !
   call ddinit(ddCosmo, n, x, y, z, rvdw)
   !
   allocate (phi(ddCosmo%ncav), psi(ddCosmo%nylm, n))
   !
   ! --------------------------   modify here  --------------------------
   !
   ! place here your favorite routine to assemble the solute's electrostatic potential
   ! and the "psi" vector. Such a routine should replace "mkrhs".
   ! for classical solutes, assembling the psi vector is straightforward; for qm solutes
   ! it requires a numerical integration (similar to the one used to compute the xc
   ! contributions in dft), as detaild in J. Chem. Phys. 141, 184108
   ! here, we compute the potential and the psi vector using the supplied routine mkrhs,
   ! which needs to be replaced by your routine.
   !
   call mkrhs(n, charge, x, y, z, ddCosmo%ncav, ddCosmo%ccav, phi, ddCosmo%nylm, psi)
   !
   ! --------------------------   end modify   --------------------------
   !
   ! now, call the ddcosmo solver
   !
   allocate (sigma(ddCosmo%nylm, n))
   !
   call cosmo(ddCosmo, .false., .true., phi, xx, psi, sigma, esolv)
   !
   if (ddCosmo%iprint.ge.3) call prtsph(ddCosmo, 'solution to the ddCOSMO equation', ddCosmo%nsph, 0, sigma)
   !
   write (6, '(1x, a, f14.6)') 'ddcosmo electrostatic solvation energy (kcal/mol):', esolv*tokcal
   !
   ! this is all for the energy. if the forces are also required, call the solver for
   ! the adjoint problem.
   ! the solution to the adjoint system is required also to compute the Fock matrix
   ! contributions.
   !
   if (ddCosmo%igrad.eq.1) then
      write(6, *)
      allocate (s(ddCosmo%nylm, n))
      allocate (fx(3, n))
      call cosmo(ddCosmo, .true., .false., xx, xx, psi, s, esolv)
      !
      if (ddCosmo%iprint.ge.3) call prtsph(ddCosmo, 'solution to the ddCOSMO adjoint equation', ddCosmo%nsph, 0, s)
      !
      ! now call the routine that computes the ddcosmo specific contributions to the forces.
      !
      call forces_dd(ddCosmo, n, phi, sigma, s, fx)
      !
      ! form the "zeta" intermediate
      !
      allocate (zeta(ddCosmo%ncav))
      call ddmkzeta(ddCosmo, s, zeta)
      !
      if (ddCosmo%iprint.ge.4) call ptcart(ddCosmo, 'zeta', ddCosmo%nsph, 0, zeta)
      !
      ! --------------------------   modify here  --------------------------
      !
      ! finally, add the contributions that depend on the derivatives of the potential
      ! on the derivatives of the potential, i.e., the electric field
      ! produced by the solute at the cavity points times the ddcosmo intermediate zeta
      ! and the solute's potential derivatives at the cavity points times zeta.
      ! the two terms are described in JCP, 141, 184108, eqs. 47, 48.
      !
      ! for a solute represented by point charges, the two terms can be rearranged as zeta
      ! contracted with the electric field of the solute plus the electric field produced
      ! by zeta, which have the physical dimension of a charge, at the nuclei times the
      ! charges. This rearrangement allows one to use fast summations methods, such as the
      ! fast multipole method, for this task.
      !
      ! a routine for point charges that follows the aformentioned strategy is provided as
      ! an example (efld). note that the coordinates should be packed into an array of
      ! dimension (3, n). we use csph, as it contains exactly this.
      ! the user will need to replace efld with his/her favorite routine.
      !
      allocate(ef(3, max(n, ddCosmo%ncav)))
      !
      ! 1. solute's electric field at the cav points times zeta:
      !
      !   compute the electric field
      !
      call efld(n, charge, ddCosmo%csph, ddCosmo%ncav, ddCosmo%ccav, ef)
      !
      !   contract it with the zeta intermediate
      !
      ii = 0
      do isph = 1, ddCosmo%nsph
         do ig = 1, ddCosmo%ngrid
            if (ddCosmo%ui(ig, isph).gt.zero) then
               ii = ii + 1
               fx(:, isph) = fx(:, isph) - zeta(ii)*ef(:, ii)
            end if
         end do
      end do
      !
      ! 2. "zeta's" electric field at the nuclei times the charges.
      !
      !   compute the "electric field"
      !
      call efld(ddCosmo%ncav, zeta, ddCosmo%ccav, n, ddCosmo%csph, ef)
      !
      !   contract it with the solute's charges.
      !
      do isph = 1, ddCosmo%nsph
         fx(:, isph) = fx(:, isph) - ef(:, isph)*charge(isph)
      end do
      !
      ! for point charges, there is no contribution from the derivatives of the psi vector.
      ! for quantum mechanical solutes, such a contribution needs to be handled via a numerical
      ! integration.
      !
      ! --------------------------   end modify   --------------------------
      !
      deallocate (zeta, ef)
      !
      if (ddCosmo%iprint.ge.1) then
         write(iout, 2000)
         2000 format(1x, 'ddCOSMO forces (atomic units):', /, &
            1x, ' atom', 15x, 'x', 15x, 'y', 15x, 'z')
         do isph = 1, ddCosmo%nsph
            write(6, '(1x, i5, 3f16.8)') isph, fx(:, isph)
         end do
      end if
   end if
   !
   ! clean up:
   !
   deallocate (x, y, z, rvdw, charge, phi, psi, sigma)
   !
   if (ddCosmo%igrad.eq.1) deallocate (s, fx)
   call memfree(ddCosmo)
   !
contains

!> Obtain the command line argument at a given index
subroutine get_argument(idx, arg)

   !> Index of command line argument, range [0:command_argument_count()]
   integer, intent(in) :: idx

   !> Command line argument
   character(len=:), allocatable, intent(out) :: arg

   integer :: length, stat

   call get_command_argument(idx, length=length, status=stat)
   if (stat /= 0) then
      return
   endif

   allocate(character(len=length) :: arg, stat=stat)
   if (stat /= 0) then
      return
   endif

   if (length > 0) then
      call get_command_argument(idx, arg, status=stat)
      if (stat /= 0) then
         deallocate(arg)
         return
      end if
   end if

end subroutine get_argument

end program main
