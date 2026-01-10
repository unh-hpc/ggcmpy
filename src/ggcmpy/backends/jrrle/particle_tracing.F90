module particle_tracing_m
   implicit none
   private

   type, public :: fields_t
      real, dimension(:,:,:), allocatable :: bx, by, bz, ex, ey, ez
      real, dimension(:), allocatable :: x_cc, y_cc, z_cc
   contains
      procedure :: init => fields_t_init
      procedure :: at => fields_t_at
   end type fields_t

contains

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

end module particle_tracing_m

module particle_tracing_f2py
   use particle_tracing_m, only: fields_t
   implicit none
   private

   type(fields_t) :: fields

   public :: load, at

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

end module particle_tracing_f2py
