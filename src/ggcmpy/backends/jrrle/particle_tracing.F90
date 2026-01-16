module particle_tracing_m
   implicit none
   private

   type, public :: fields_t
      real, dimension(:,:,:), allocatable :: bx, by, bz, ex, ey, ez
      real, dimension(:), allocatable :: x_cc, y_cc, z_cc
      real, dimension(:), allocatable :: x_nc, y_nc, z_nc
   contains
      procedure :: init => fields_t_init
      procedure :: at => fields_t_at
      procedure :: interpolate => fields_t_interpolate
      procedure :: interpolate_yee => fields_t_interpolate_yee
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

   subroutine fields_t_init(this, bx, by, bz, ex, ey, ez, x_cc, y_cc, z_cc, x_nc, y_nc, z_nc)
      class(fields_t), intent(inout) :: this
      real, dimension(:,:,:), intent(in) :: bx, by, bz, ex, ey, ez
      real, dimension(:), intent(in) :: x_cc, y_cc, z_cc
      real, dimension(:), intent(in), optional :: x_nc, y_nc, z_nc
      this%bx = bx
      this%by = by
      this%bz = bz
      this%ex = ex
      this%ey = ey
      this%ez = ez
      this%x_cc = x_cc
      this%y_cc = y_cc
      this%z_cc = z_cc
      if (present(x_nc)) then
         this%x_nc = x_nc
         this%y_nc = y_nc
         this%z_nc = z_nc
      end if
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

   real function ip(x, y, z, crd_x, crd_y, crd_z, field)
      real, intent(in) :: x, y, z
      real, dimension(:), intent(in) :: crd_x, crd_y, crd_z
      real, dimension(:,:,:), intent(in) :: field
      integer :: i, j, k
      real :: xd, yd, zd
      real :: c000, c100, c010, c110, c001, c101, c011, c111
      real :: x0, x1, y0, y1, z0, z1

      ! Find indices
      i = find_index(crd_x, x)
      j = find_index(crd_y, y)
      k = find_index(crd_z, z)

      ! Get cell bounds
      x0 = crd_x(i)
      x1 = crd_x(i+1)
      y0 = crd_y(j)
      y1 = crd_y(j+1)
      z0 = crd_z(k)
      z1 = crd_z(k+1)

      ! Compute normalized distances
      xd = (x - x0) / (x1 - x0)
      yd = (y - y0) / (y1 - y0)
      zd = (z - z0) / (z1 - z0)

      ! Get field values at corners
      c000 = field(i  , j  , k  )
      c100 = field(i+1, j  , k  )
      c010 = field(i  , j+1, k  )
      c110 = field(i+1, j+1, k  )
      c001 = field(i  , j  , k+1)
      c101 = field(i+1, j  , k+1)
      c011 = field(i  , j+1, k+1)
      c111 = field(i+1, j+1, k+1)

      ! Trilinear interpolation
      ip = &
         c000 * (1-xd)*(1-yd)*(1-zd) + &
         c100 * xd    *(1-yd)*(1-zd) + &
         c010 * (1-xd)*yd    *(1-zd) + &
         c110 * xd    *yd    *(1-zd) + &
         c001 * (1-xd)*(1-yd)*zd     + &
         c101 * xd    *(1-yd)*zd     + &
         c011 * (1-xd)*yd    *zd     + &
         c111 * xd    *yd    *zd
   end function ip

   real function fields_t_interpolate(this, x, y, z, m)
      class(fields_t), target, intent(in) :: this
      real, intent(in) :: x, y, z
      integer, intent(in) :: m

      real, dimension(:,:,:), pointer :: field
      integer :: i, j, k
      real :: xd, yd, zd
      real :: c000, c100, c010, c110, c001, c101, c011, c111
      real :: x0, x1, y0, y1, z0, z1

      select case (m)
       case (0)
         field => this%bx
       case (1)
         field => this%by
       case (2)
         field => this%bz
       case (3)
         field => this%ex
       case (4)
         field => this%ey
       case (5)
         field => this%ez
       case default
         stop 'Invalid field index in fields_t_interpolate'
      end select

      fields_t_interpolate = ip(x, y, z, this%x_cc, this%y_cc, this%z_cc, field)
   end function fields_t_interpolate

   real function fields_t_interpolate_yee(this, x, y, z, m)
      class(fields_t), intent(in) :: this
      real, intent(in) :: x, y, z
      integer, intent(in) :: m

      select case(m)
       case (0)
         fields_t_interpolate_yee = ip(x, y, z, this%x_nc, this%y_cc, this%z_cc, this%bx)
       case (1)
         fields_t_interpolate_yee = ip(x, y, z, this%x_cc, this%y_nc, this%z_cc, this%by)
       case (2)
         fields_t_interpolate_yee = ip(x, y, z, this%x_cc, this%y_cc, this%z_nc, this%bz)
       case (3)
         fields_t_interpolate_yee = ip(x, y, z, this%x_cc, this%y_nc, this%z_nc, this%ex)
       case (4)
         fields_t_interpolate_yee = ip(x, y, z, this%x_nc, this%y_cc, this%z_nc, this%ey)
       case (5)
         fields_t_interpolate_yee = ip(x, y, z, this%x_nc, this%y_nc, this%z_cc, this%ez)
       case default
         stop 'Invalid field index in fields_t_interpolate_yee'
      end select
   end function fields_t_interpolate_yee

   subroutine boris_integrator_t_init(this, q, m)
      class(boris_integrator_t), intent(inout) :: this
      real, intent(in) :: q  ! charge [C]
      real, intent(in) :: m  ! mass [kg]
      this%q = q
      this%m = m
   end subroutine boris_integrator_t_init

   subroutine boris_integrator_t_integrate(this, x0, v0, get_E, get_B, t_max, dt_max, gyro_max, data, n_out)
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
      real, intent(in) :: t_max, dt_max, gyro_max
      real, dimension(:, 0:), intent(out) :: data
      integer, intent(out) :: n_out

      integer :: step, n_data
      real :: t
      real, dimension(3) :: x, v, E, B
      real :: qprime, om_c, dt
      real, dimension(3) :: h, s
      real, parameter :: pi = 3.14159265358979323846

      n_data = size(data, 2)

      t = 0.0
      x = x0
      v = v0
      qprime = 0.5 * this%q / this%m
      ! times, positions, velocities = [], [], []
      B = get_B(x)
      step = 0
      do while (t < t_max)
         if (step < n_data) then
            data(1, step) = t
            data(2:4, step) = x
            data(5:7, step) = v
            step = step + 1
         end if

         om_c = norm2(2. * qprime * B)  ! gyro frequency
         dt = min(dt_max, gyro_max * 2.0 * pi / om_c)

         x = x + 0.5 * dt * v
         B = get_B(x)
         E = get_E(x)
         v = v + dt * qprime * E
         h = dt * qprime * B
         s = 2. * h / (1. + norm2(h) ** 2)
         v = v + cross(v + cross(v, h), s)
         v = v + dt * qprime * E
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
   public :: load_yee, interpolate_yee
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

   subroutine load_yee(nx, ny, nz, nx_nc, ny_nc, nz_nc, &
      bx, by, bz, ex, ey, ez, &
      x_cc, y_cc, z_cc, x_nc, y_nc, z_nc)
      integer, intent(in) :: nx, ny, nz, nx_nc, ny_nc, nz_nc
      !f2py intent(hide) :: nx, ny, nz, nx_nc, ny_nc, nz_nc
      real, dimension(nx_nc,ny,nz), intent(in) :: bx
      real, dimension(nx,ny_nc,nz), intent(in) :: by
      real, dimension(nx,ny,nz_nc), intent(in) :: bz
      real, dimension(nx,ny_nc,nz_nc), intent(in) :: ex
      real, dimension(nx_nc,ny,nz_nc), intent(in) :: ey
      real, dimension(nx_nc,ny_nc,nz), intent(in) :: ez
      real, dimension(nx), intent(in) :: x_cc
      real, dimension(ny), intent(in) :: y_cc
      real, dimension(nz), intent(in) :: z_cc
      real, dimension(nx_nc), intent(in) :: x_nc
      real, dimension(ny_nc), intent(in) :: y_nc
      real, dimension(nz_nc), intent(in) :: z_nc

      call fields%init(bx, by, bz, ex, ey, ez, x_cc, y_cc, z_cc, x_nc, y_nc, z_nc)
   end subroutine load_yee

   real function at(i, j, k, m)
      integer, intent(in) :: i, j, k, m

      at = fields%at(i, j, k, m)
   end function at

   real function interpolate(x, y, z, m)
      real, intent(in) :: x, y, z
      integer, intent(in) :: m

      interpolate = fields%interpolate(x, y, z, m)
   end function interpolate

   real function interpolate_yee(x, y, z, m)
      real, intent(in) :: x, y, z
      integer, intent(in) :: m

      interpolate_yee = fields%interpolate_yee(x, y, z, m)
   end function interpolate_yee


   subroutine boris_init(q, m)
      real, intent(in) :: q  ! charge [C]
      real, intent(in) :: m  ! mass [kg]

      call boris_integrator%init(q, m)
   end subroutine boris_init

   function get_B(x) result(B)
      real, dimension(3), intent(in) :: x
      real, dimension(3) :: B
      B = [interpolate_yee(x(1), x(2), x(3), 0), interpolate_yee(x(1), x(2), x(3), 1), interpolate_yee(x(1), x(2), x(3), 2)]
   end function get_B

   function get_E(x) result(E)
      real, dimension(3), intent(in) :: x
      real, dimension(3) :: E
      E = [interpolate_yee(x(1), x(2), x(3), 3), interpolate_yee(x(1), x(2), x(3), 4), interpolate_yee(x(1), x(2), x(3), 5)]
   end function get_E

   subroutine boris_integrate(x0, v0, t_max, dt_max, gyro_max, data, n_out, n_data)
      real, dimension(3), intent(in) :: x0
      real, dimension(3), intent(in) :: v0
      real, intent(in) :: t_max, dt_max, gyro_max
      real, dimension(7, n_data) :: data
      integer, intent(out) :: n_out
      integer, intent(in) :: n_data
      !f2py intent(hide) :: n_data

      call boris_integrator%integrate(x0, v0, get_E, get_B, t_max, dt_max, gyro_max, data, n_out)
   end subroutine boris_integrate

end module particle_tracing_f2py
