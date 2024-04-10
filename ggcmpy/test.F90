module module_calculate_stats
   use iso_fortran_env
   use iso_c_binding
   implicit none

contains
   subroutine calculate_sum(array, array_sum, n) bind(C)
      integer(c_int), value :: n
      real(kind=c_float), intent(in), dimension(n)::array
      real(kind=c_float), intent(out)::array_sum
      array_sum=sum(array)
   end subroutine calculate_sum

   subroutine calculate_mean(array, array_mean, n) bind(C)
      integer(c_int), value :: n
      real(kind=c_float), intent(in), dimension(n)::array
      real(kind=c_float), intent(out)::array_mean
      array_mean=sum(array)/n
   end subroutine calculate_mean

end module module_calculate_stats
