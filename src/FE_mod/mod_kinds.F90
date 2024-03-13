module mod_kinds

   use iso_fortran_env

   implicit none

   private
   public :: wp, zero, one_half, one, two, three, four

#ifdef REAL32
   integer, parameter :: wp = real32   !! default real kind [4 bytes]
#elif REAL64
   integer, parameter :: wp = real64   !! default real kind [8 bytes]
#elif REAL128
   integer, parameter :: wp = real128  !! default real kind [16 bytes]
#else
   integer, parameter :: wp = real64   !! default real kind [8 bytes]
#endif

   ! parameters:
   real(wp), parameter :: zero = 0.0_wp
   real(wp), parameter :: one_half = 0.5_wp
   real(wp), parameter :: one = 1.0_wp
   real(wp), parameter :: two = 2.0_wp
   real(wp), parameter :: three = 3.0_wp
   real(wp), parameter :: four = 4.0_wp
   real(wp), parameter :: pi = acos(-one)

end module mod_kinds
