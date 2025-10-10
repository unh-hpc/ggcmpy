
module f2py
   ! wrappers for OpenGGCM routines that can not be directly wrapped by f2py
   ! due to use of Fortran types, etc
   use iono_psol_m, only: iono_potential_solve_t
   implicit none
   private

   public :: epoch1966, cotr
   public :: iono_potential_solve_initialize, iono_potential_solve_finalize
   public :: iono_potential_solve_setup, iono_potential_solve

   type(iono_potential_solve_t), save :: global_solver
   logical, save :: global_solver_initialized = .false.

contains

   subroutine epoch1966(dsecs, iy, mo, id, ih, mi, sec)
      use date_lib, only: epoch1966_ => epoch1966
      real*8, intent(out) :: dsecs
      !f2py intent(out) :: dsecs
      integer, intent(in) :: iy, mo, id, ih, mi
      real, intent(in) :: sec

      call epoch1966_(dsecs, iy, mo, id, ih, mi, sec, 0)
   end subroutine epoch1966

   subroutine cotr(dsecs, cfr, cto, r1, r2)
      use new_cotr, only: cotr_set, cotr_ => cotr, transform
      real*8, intent(in) :: dsecs
      character*3, intent(in) :: cfr, cto
      real, dimension(3), intent(in) :: r1
      real, dimension(3), intent(out) :: r2

      type(transform) :: tran

      call cotr_set(dsecs, tran)
      call cotr_(tran, cfr, cto, r1(1), r1(2), r1(3), r2(1), r2(2), r2(3))
   end subroutine cotr

   subroutine iono_potential_solve_initialize(nphi, nthe, nx, ny, ngord, npany)
      use iono_psol_m, only: iono_potential_solve_initialize_ => iono_potential_solve_initialize
      integer, intent(in) :: nphi, nthe, nx, ny, ngord, npany

      if (global_solver_initialized) then
         call iono_potential_solve_finalize()
      end if
      call iono_potential_solve_initialize_(global_solver, nphi, nthe, nx, ny, ngord, npany)
      global_solver_initialized = .true.
   end subroutine iono_potential_solve_initialize

   subroutine iono_potential_solve_finalize()
      use iono_psol_m, only: iono_potential_solve_finalize_ => iono_potential_solve_finalize

      if (global_solver_initialized) then
         call iono_potential_solve_finalize_(global_solver)
         global_solver_initialized = .false.
      end if
   end subroutine iono_potential_solve_finalize

   subroutine iono_potential_solve_setup(np, nt, sigp, sigh)
      use iono_psol_m, only: iono_potential_solve_setup_ => iono_potential_solve_setup
      use constants_m, only: sig0

      integer :: np, nt
      !f2py intent(hide) :: np, nt
      real, dimension(np, nt), intent(in) :: sigp, sigh

      real, dimension(np, nt) :: sigp_scaled, sigh_scaled

      sigp_scaled = sigp / sig0
      sigh_scaled = sigh / sig0

      if (.not. global_solver_initialized) then
         write(*,*) 'Error: global iono_potential_solve not initialized.'
         stop
      end if

      call iono_potential_solve_setup_(global_solver, sigp_scaled, sigh_scaled)
   end subroutine iono_potential_solve_setup

   subroutine iono_potential_solve(np, nt, xjpa, pot)
      use iono_psol_m, only: iono_potential_solve_ => iono_potential_solve
      use constants_m, only: phi0, j0

      integer :: np, nt
      !f2py intent(hide) :: np, nt
      real, dimension(np, nt), intent(in) :: xjpa
      real, dimension(np, nt), intent(out) :: pot

      real, dimension(np, nt) :: xjpa_scaled

      if (.not. global_solver_initialized) then
         write(*,*) 'Error: global iono_potential_solve not initialized.'
         stop
      end if

      xjpa_scaled = xjpa / j0
      call iono_potential_solve_(global_solver, xjpa_scaled, pot)

      pot = pot * phi0
   end subroutine iono_potential_solve

end module f2py
