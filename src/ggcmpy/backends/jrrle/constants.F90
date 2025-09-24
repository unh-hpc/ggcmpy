module constants_m
   implicit none
   private

   real, parameter, public :: pi = 3.1415927, pi2 = 6.2831855

   real, parameter, public :: sig0 = 0.11889699
   real, parameter, public :: phi0= 0.13037126E+10
   real, parameter, public :: j0 = 0.38188455E-05
   real, parameter, public :: n0 = 0.10000000E+11
   real, parameter, public :: p0 = 0.37193298E-03
   real, parameter, public :: v0 = 6692976.5

   ! MPI tags
   integer, parameter, public :: ID_IONEGO = 1000
   integer, parameter, public :: ID_IOALL = 2000
   integer, parameter, public :: ID_POT = 3000

   integer, parameter, public :: tty = 0

end module constants_m
