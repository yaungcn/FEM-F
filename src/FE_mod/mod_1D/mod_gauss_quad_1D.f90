module mod_gauss_quad_1D
   use quadrature_module, wp => quadrature_wp
   use mod_basis_func_1D
   use mod_co_func_1D
   implicit none

   private
   public :: quad_1D, quad_A_1D, quad_b_1D

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
      !! TODO: used in assemble matrix A.
      !! SUIT: 1D linear. Dirchlet boundary condition.
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
      !! TODO: used in assemble vector b.
      !! SUIT: 1D linear. Dirchlet boundary condition.
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

   function quad_solution_1D(self, x) result(f)
      !! TODO: used to quadrature the solution function.
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

   end function quad_solution_1D

end module mod_gauss_quad_1D
