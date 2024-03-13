module mod_co_func_2D
   use mod_kinds
   implicit none

   private
   public :: function_c, function_f, function_g

   abstract interface
   !! add pure
      pure function cofunc_1d(x)
         import :: wp
         real(wp), intent(in) :: x
         real(wp) :: cofunc_1d
      end function cofunc_1d
   !! add pure, 2D
      pure function cofunc_2d(x, y)
         import :: wp
         real(wp), intent(in) :: x, y
         real(wp) :: cofunc_2d
      end function cofunc_2d
   !! add pure, 2D, time
      pure function cofunc_2d_t(x, y, t)
         import :: wp
         real(wp), intent(in) :: x, y, t
         real(wp) :: cofunc_2d_t
      end function cofunc_2d_t
   end interface

contains
!! add pure
   pure function function_c(x, y)
      real(wp), intent(in) :: x, y
      real(wp) :: function_c
      function_c = exp(x)
   end function function_c

   pure function function_f(x, y)
      real(wp), intent(in) :: x, y
      real(wp) :: function_f
      function_f = -exp(x)*(cos(x) - 2*sin(x) - x*cos(x) - x*sin(x))
   end function function_f

   pure function function_g(x, y)
      real(wp), intent(in) :: x, y
      real(wp) :: function_g
      select case (int(x))
      case (0)
         function_g = 0
      case (1)
         function_g = cos(1.0_wp)
      end select
   end function function_g
end module mod_co_func_2D
