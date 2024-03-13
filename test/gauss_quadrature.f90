!************************************************************************************
!>
!  Test of the integration routine.
!  Integrates various functions and prints results.

program main
   use mod_kinds
   use quadrature_module

   implicit none

   real(wp), parameter :: zero = 0.0_wp
   real(wp), parameter :: one = 1.0_wp
   real(wp), parameter :: two = 2.0_wp
   real(wp), parameter :: three = 3.0_wp
   real(wp), parameter :: pi = acos(-one)
   real(wp), parameter :: tol = 100000*epsilon(one) !! error tolerance
   real(wp), parameter :: abs_error_tol_for_check = 100*tol !! tol for pass/fail check

   type, extends(integration_class_1d) :: sin_type
      integer  :: ifunc = 0     !! which function to use
      real(wp) :: amp = zero  !! amplitude
      real(wp) :: freq = zero  !! frequency
      real(wp) :: phase = zero  !! phase
      integer  :: n_evals = 0     !! number of function evaluations
   end type sin_type

   type, extends(integration_class_2d) :: my_doub
      integer :: ifunc = 0   !! which function to use
      integer :: n_evals = 0   !! number of function evaluations
   end type my_doub
   type, extends(integration_class_3d) :: my_trip
      integer :: ifunc = 0   !! which function to use
      integer :: n_evals = 0   !! number of function evaluations
   end type my_trip

   type(my_doub)  :: doub

   real(wp) :: ans, err, answer
   real(wp) :: xl, xu, yl, yu
   integer  :: ierr    !! error code
   integer  :: i       !! counter
   integer  :: meth    !! method number
   integer  :: itest   !! test number for output printing

   itest = 0

   ! csv header:
   write (*, '(A)') 'Test number, Dimension, Method, tol, ans, ierr, err, Function evaluations, Actual Error'

   do i = 1, size(set_of_quadrature_methods)  !test all the methods

      meth = set_of_quadrature_methods(i)%n_points

      !============================================
      ! double integral tests
      !============================================

      doub%ifunc = 1
      xl = 0.3827812_wp
      xu = 1.02928_wp
      yl = 22.91281_wp
      yu = -111.928_wp
      call run_2d_test()

      doub%ifunc = 2
      xl = 0.0_wp
      xu = pi
      yl = 0.0_wp
      yu = 2.0_wp*pi
      call run_2d_test()

   end do

contains
!************************************************************************************

   !*************************************************************
   subroutine run_2d_test()

      implicit none

      itest = doub%ifunc

      !set up the class
      call doub%initialize(fxy=test_2d_func, xl=xl, xu=xu, yl=yl, &
         yu=yu, tolx=tol, toly=tol, methodx=meth, methody=meth)

      !reset number of function evaluations:
      doub%n_evals = 0

      !integrate the function:
      call doub%integrate(ans, ierr, err)

      !get the true answer:
      answer = test_2d_integral(doub, xl, xu, yl, yu)

      !print results:
      write (*, '(1p,I3,A,A,A,A,A,E15.5,A,E15.5,A,I5,A,E15.5,A,I11,A,E15.5)') &
         itest, ',', &
         '2D', ',', &
         trim(set_of_quadrature_methods(i)%name), ',', &
         tol, ',', &
         ans, ',', &
         ierr, ',', &
         err, ',', &
         doub%n_evals, ',', &
         answer - ans

      if (abs(answer - ans) > abs_error_tol_for_check) error stop 'TEST FAILED'

   end subroutine run_2d_test
   !*************************************************************

   !*************************************************************
   function test_2d_func(me, x, y) result(f)

      !! The function is f(x,y)

      implicit none

      class(integration_class_2d), intent(inout)   :: me
      real(wp), intent(in)  :: x
      real(wp), intent(in)  :: y
      real(wp)              :: f

      select type (me)
       class is (my_doub)

         select case (me%ifunc)
          case (1)

            f = x*y

          case (2)

            f = sin(x) + cos(y)

          case default
            error stop 'Error in test_2d_func: invalid value of ifunc'
         end select

         me%n_evals = me%n_evals + 1

       class default
         error stop 'Error in test_2d_func: invalid class.'
      end select

   end function test_2d_func
   !*************************************************************

   !*************************************************************
   function test_2d_integral(me, xl, xu, yl, yu) result(f)

      !! The double integral of f(x,y) dx dy = x*y from xl->xu, yl->yu

      implicit none

      class(integration_class_2d), intent(inout)  :: me
      real(wp), intent(in)  :: xl
      real(wp), intent(in)  :: xu
      real(wp), intent(in)  :: yl
      real(wp), intent(in)  :: yu
      real(wp)              :: f

      select type (me)
       class is (my_doub)

         select case (me%ifunc)
          case (1)

            f = (xu**2/two - xl**2/two)*(yu**2/two - yl**2/two)

          case (2)

            f = (-cos(xu) + cos(xl))*(yu - yl) + (xu - xl)*(sin(yu) - sin(yl))

          case default
            error stop 'Error in test_2d_integral: invalid value of ifunc'
         end select

       class default
         error stop 'Error in test_2d_integral: invalid class.'
      end select

   end function test_2d_integral
   !*************************************************************

!************************************************************************************
end program main
!************************************************************************************
