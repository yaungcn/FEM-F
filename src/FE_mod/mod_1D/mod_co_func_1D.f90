module mod_co_func_1D
   use quadrature_module, only: wp => quadrature_wp
   implicit none

   private
   public :: func, func_a, func_f, func_g

   type, abstract :: func
   contains
      procedure(func_abstract), nopass, deferred :: func_pass
   end type

   abstract interface
      function func_abstract(x)
         import :: wp
         real(wp), intent(in) :: x
         real(wp) :: func_abstract
      end function func_abstract
   end interface

   type, extends(func) :: func_a
   contains
      procedure, nopass :: func_pass => function_a
   end type

   type, extends(func) :: func_f
   contains
      procedure, nopass :: func_pass => function_f
   end type

   type, extends(func) :: func_g
   contains
      procedure, nopass :: func_pass => function_g
   end type
contains

   function function_a(x)
      real(wp), intent(in) :: x
      real(wp) :: function_a
      function_a = exp(x)
   end function function_a

   function function_f(x)
      real(wp), intent(in) :: x
      real(wp) :: function_f
      function_f = -exp(x)*(cos(x) - 2*sin(x) - x*cos(x) - x*sin(x))
   end function function_f

   function function_g(x)
      real(wp), intent(in) :: x
      real(wp) :: function_g
      select case (int(x))
      case (0)
         function_g = 0
      case (1)
         function_g = cos(1.0_wp)
      end select
   end function function_g
end module mod_co_func_1D
