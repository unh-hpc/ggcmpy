module particle_tracing_m
   implicit none
   private

   type, public :: fields_t
      real, dimension(:,:,:), allocatable :: bx, by, bz, ex, ey, ez
      real, dimension(:), allocatable :: x_cc, y_cc, z_cc
   contains
      procedure :: init => fields_t_init
      procedure :: at => fields_t_at
      procedure :: interpolate => fields_t_interpolate
   end type fields_t

   type, public :: boris_integrator_t
      real :: q  ! charge [C]
      real :: m  ! mass [kg]
   contains
      procedure :: init => boris_integrator_t_init
      procedure :: integrate => boris_integrator_t_integrate
   end type boris_integrator_t

contains

   integer function find_index(crd, val)
      real, dimension(:), intent(in) :: crd
      real, intent(in) :: val
      integer :: n, left, right, mid

      n = size(crd)
      if (val < crd(1)) then
         find_index = 1
         return
      else if (val >= crd(n)) then
         find_index = n - 1
         return
      end if

      left = 1
      right = n - 1
      do while (left <= right)
         mid = (left + right) / 2
         if (val < crd(mid)) then
            right = mid - 1
         else if (val >= crd(mid+1)) then
            left = mid + 1
         else
            find_index = mid
            return
         end if
      end do
      stop 'should not get here'     ! Should not reach here if input is valid
   end function find_index

   subroutine fields_t_init(this, bx, by, bz, ex, ey, ez, x_cc, y_cc, z_cc)
      class(fields_t), intent(inout) :: this
      real, dimension(:,:,:), intent(in) :: bx, by, bz, ex, ey, ez
      real, dimension(:), intent(in) :: x_cc, y_cc, z_cc
      this%bx = bx
      this%by = by
      this%bz = bz
      this%ex = ex
      this%ey = ey
      this%ez = ez
      this%x_cc = x_cc
      this%y_cc = y_cc
      this%z_cc = z_cc
   end subroutine fields_t_init

   real function fields_t_at(this, i, j, k, m)
      class(fields_t), intent(in) :: this
      integer, intent(in) :: i, j, k, m

      select case (m)
       case (0)
         fields_t_at = this%bx(i,j,k)
       case (1)
         fields_t_at = this%by(i,j,k)
       case (2)
         fields_t_at = this%bz(i,j,k)
       case (3)
         fields_t_at = this%ex(i,j,k)
       case (4)
         fields_t_at = this%ey(i,j,k)
       case (5)
         fields_t_at = this%ez(i,j,k)
       case default
         stop 'Invalid field index in fields_t_at'
      end select
   end function fields_t_at

   real function fields_t_interpolate(this, x, y, z, m)
      class(fields_t), intent(in) :: this
      real, intent(in) :: x, y, z
      integer, intent(in) :: m
      integer :: i, j, k
      real :: xd, yd, zd
      real :: c000, c100, c010, c110, c001, c101, c011, c111
      real :: x0, x1, y0, y1, z0, z1

      ! Find indices
      i = find_index(this%x_cc, x)
      j = find_index(this%y_cc, y)
      k = find_index(this%z_cc, z)

      ! Get cell bounds
      x0 = this%x_cc(i)
      x1 = this%x_cc(i+1)
      y0 = this%y_cc(j)
      y1 = this%y_cc(j+1)
      z0 = this%z_cc(k)
      z1 = this%z_cc(k+1)

      ! Compute normalized distances
      xd = (x - x0) / (x1 - x0)
      yd = (y - y0) / (y1 - y0)
      zd = (z - z0) / (z1 - z0)

      ! Get field values at corners
      c000 = this%at(i  , j  , k  , m)
      c100 = this%at(i+1, j  , k  , m)
      c010 = this%at(i  , j+1, k  , m)
      c110 = this%at(i+1, j+1, k  , m)
      c001 = this%at(i  , j  , k+1, m)
      c101 = this%at(i+1, j  , k+1, m)
      c011 = this%at(i  , j+1, k+1, m)
      c111 = this%at(i+1, j+1, k+1, m)

      ! Trilinear interpolation
      fields_t_interpolate = &
         c000 * (1-xd)*(1-yd)*(1-zd) + &
         c100 * xd    *(1-yd)*(1-zd) + &
         c010 * (1-xd)*yd    *(1-zd) + &
         c110 * xd    *yd    *(1-zd) + &
         c001 * (1-xd)*(1-yd)*zd     + &
         c101 * xd    *(1-yd)*zd     + &
         c011 * (1-xd)*yd    *zd     + &
         c111 * xd    *yd    *zd

   end function fields_t_interpolate

   subroutine boris_integrator_t_init(this, q, m)
      class(boris_integrator_t), intent(inout) :: this
      real, intent(in) :: q  ! charge [C]
      real, intent(in) :: m  ! mass [kg]
      this%q = q
      this%m = m
   end subroutine boris_integrator_t_init

   subroutine boris_integrator_t_integrate(this, x0, v0, get_E, get_B, t_max, dt, data, n_out)
      class(boris_integrator_t), intent(in) :: this
      real, dimension(3), intent(in) :: x0
      real, dimension(3), intent(in) :: v0
      interface
         function get_E(x) result(E)
            real, dimension(3), intent(in) :: x
            real, dimension(3) :: E
         end function get_E
         function get_B(x) result(B)
            real, dimension(3), intent(in) :: x
            real, dimension(3) :: B
         end function get_B
      end interface
      real, intent(in) :: t_max
      real, intent(in) :: dt
      real, dimension(:, 0:), intent(out) :: data
      integer, intent(out) :: n_out

      integer :: step, n_data
      real :: t
      real, dimension(3) :: x, v, E, B
      real :: qprime
      real, dimension(3) :: h, s

      n_data = size(data, 2)

      t = 0.0
      x = x0
      v = v0
      qprime = 0.5 * dt * this%q / this%m
      ! times, positions, velocities = [], [], []
      step = 0
      do while (t < t_max)
         if (step < n_data) then
            data(1, step) = t
            data(2:4, step) = x
            data(5:7, step) = v
            step = step + 1
         end if
         ! times.append(t)
         ! positions.append(x.copy())
         ! velocities.append(v.copy())
         B = get_B(x)
         E = get_E(x)
         x = x + 0.5 * dt * v
         v = v + qprime * E
         h = qprime * B
         s = 2. * h / (1. + abs(h) ** 2)
         v = v + cross(v + cross(v, h), s)
         v = v + qprime * E
         x = x + 0.5 * dt * v
         t = t + dt
      end do

      if (step < n_data) then
         data(1, step) = t
         data(2:4, step) = x
         data(5:7, step) = v
         step = step + 1
      end if
      n_out = step

   contains
      function cross(a, b) result(c)
         real, dimension(3), intent(in) :: a, b
         real, dimension(3) :: c
         c(1) = a(2)*b(3) - a(3)*b(2)
         c(2) = a(3)*b(1) - a(1)*b(3)
         c(3) = a(1)*b(2) - a(2)*b(1)
      end function cross
   end subroutine boris_integrator_t_integrate

end module particle_tracing_m


module particle_tracing_f2py
   use particle_tracing_m, only: fields_t, boris_integrator_t
   implicit none
   private

   type(fields_t) :: fields
   type(boris_integrator_t) :: boris_integrator

   public :: load, at, interpolate
   public :: boris_init, boris_integrate

contains

   subroutine load(nx, ny, nz, bx, by, bz, ex, ey, ez, x_cc, y_cc, z_cc)
      integer, intent(in) :: nx, ny, nz
      !f2py intent(hide) :: nx, ny, nz
      real, dimension(nx,ny,nz), intent(in) :: bx, by, bz, ex, ey, ez
      real, dimension(nx), intent(in) :: x_cc
      real, dimension(ny), intent(in) :: y_cc
      real, dimension(nz), intent(in) :: z_cc

      call fields%init(bx, by, bz, ex, ey, ez, x_cc, y_cc, z_cc)
   end subroutine load

   real function at(i, j, k, m)
      integer, intent(in) :: i, j, k, m

      at = fields%at(i, j, k, m)
   end function at

   real function interpolate(x, y, z, m)
      real, intent(in) :: x, y, z
      integer, intent(in) :: m

      interpolate = fields%interpolate(x, y, z, m)
   end function interpolate

   subroutine boris_init(q, m)
      real, intent(in) :: q  ! charge [C]
      real, intent(in) :: m  ! mass [kg]

      call boris_integrator%init(q, m)
   end subroutine boris_init

   function get_B(x) result(B)
      real, dimension(3), intent(in) :: x
      real, dimension(3) :: B
      B = [interpolate(x(1), x(2), x(3), 0), interpolate(x(1), x(2), x(3), 1), interpolate(x(1), x(2), x(3), 2)]
   end function get_B

   function get_E(x) result(E)
      real, dimension(3), intent(in) :: x
      real, dimension(3) :: E
      E = [interpolate(x(1), x(2), x(3), 3), interpolate(x(1), x(2), x(3), 4), interpolate(x(1), x(2), x(3), 5)]
   end function get_E

   subroutine boris_integrate(x0, v0, t_max, dt, data, n_out, n_data)
      real, dimension(3), intent(in) :: x0
      real, dimension(3), intent(in) :: v0
      real, intent(in) :: t_max
      real, intent(in) :: dt
      real, dimension(7, n_data) :: data
      integer, intent(out) :: n_out
      integer, intent(in) :: n_data
      !f2py intent(hide) :: n_data

      call boris_integrator%integrate(x0, v0, get_E, get_B, t_max, dt, data, n_out)
   end subroutine boris_integrate

end module particle_tracing_f2py
