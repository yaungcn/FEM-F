program main
   use quadrature_module, wp => quadrature_wp

   implicit none

   type(integration_class_1d) :: quad_1D
   ! type(integration_class_2d) :: quad_2D
   ! type(integration_class_3d) :: quad_3D

   real(wp) :: a, b, ans, err
   real(wp), parameter :: tol = 10000*epsilon(1.0_wp)
   integer :: meth = 10 ! gauss point 6-14
   integer :: ierr

   a = 0.0_wp
   b = 4*ATAN(1.0_wp)

   call quad_1D%initialize(fx=test_func, xl=a, xu=b, tolx=tol, methodx=meth)

   call quad_1D%integrate(ans, ierr, err)

   write (*, '(A)') "ans"
   write (*, '(E15.5)') ans
contains

   function test_func(self, x) result(f)
      class(integration_class_1d), intent(inout) :: self
      real(wp), intent(in) :: x
      real(wp) :: f

      f = cos(x)

   end function test_func
end program main

