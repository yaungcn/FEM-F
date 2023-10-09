module mod_gauss_quad
   use quadrature_module, wp => quadrature_wp
   use mod_basis_func
   use mod_co_func
   implicit none

   private
   public :: quad_1D, quad_2D, quad_A_1D, quad_b_1D

   type, extends(integration_class_1d) :: quad_1D
      private
      class(func), pointer :: co_func
            !! The coefficient function of the integrand.
      class(local_basis), pointer :: basis_trial
            !! local trial basis function
      class(local_basis), pointer :: basis_test
            !! local trial basis function

   contains
      procedure :: init => init_1D
      ! procedure :: quad_func => quad_func_1D
   end type quad_1D

   type, extends(integration_class_2d) :: quad_2D

   end type quad_2D

contains
   subroutine init_1D(self, co_func, &
                      basis_trial, basis_test)
      class(quad_1D), intent(inout) :: self
      class(func), target :: co_func
      class(local_basis), target :: basis_trial, basis_test

      self%co_func => co_func
      self%basis_trial => basis_trial
      self%basis_test => basis_test

   end subroutine init_1D

   function quad_A_1D(self, x) result(f)
      class(integration_class_1d), intent(inout) :: self
      real(wp), intent(in) :: x
      real(wp) :: f

      select type (self)
      class is (quad_1D)
         f = self%co_func%func_pass(x)* &
             self%basis_trial%basis_f(x)* &
             self%basis_test%basis_f(x)
      class default
         error stop "Undefined Class"
      end select

   end function quad_A_1D

   function quad_b_1D(self, x) result(f)
      class(integration_class_1d), intent(inout) :: self
      real(wp), intent(in) :: x
      real(wp) :: f

      select type (self)
      class is (quad_1D)
         f = self%co_func%func_pass(x)* &
             self%basis_test%basis_f(x)
      class default
         error stop "Undefined Class"
      end select

   end function quad_b_1D
end module mod_gauss_quad
