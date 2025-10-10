
module f2py
   ! wrappers for OpenGGCM routines that can not be directly wrapped by f2py
   ! due to use of Fortran types, etc
   implicit none
   private

   public :: epoch1966, cotr

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

end module f2py
