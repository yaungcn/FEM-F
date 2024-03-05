module mod_gauss_quad_2D
   use quadrature_module, wp => quadrature_wp
   use mod_field_2D, only: field
   use mod_local_basis_func_2D
   implicit none

   private
   public :: func_to_quad_2d, afunc, bfunc

   abstract interface
      pure function local_basis_func(x, y, local_vertices, basis_type, basis_index, der_x, der_y) result(func)
         import :: wp
         real(wp), intent(in) :: x, y
         real(wp), intent(in) :: local_vertices(2, 3)
         integer, intent(in) :: basis_type, basis_index, der_x, der_y
         real(wp) :: func
      end function local_basis_func
      !! 1D
      pure function cofunc_1d(x)
         import :: wp
         real(wp), intent(in) :: x
         real(wp) :: cofunc_1d
      end function cofunc_1d
      !! 2D
      pure function cofunc_2d(x, y)
         import :: wp
         real(wp), intent(in) :: x, y
         real(wp) :: cofunc_2d
      end function cofunc_2d
      !! 2D, time
      pure function cofunc_2d_t(x, y, t)
         import :: wp
         real(wp), intent(in) :: x, y, t
         real(wp) :: cofunc_2d_t
      end function cofunc_2d_t
   end interface

   type, extends(integration_class_2d), public :: func_to_quad_2d

      private
      real(wp) :: vertices(2, 3)
      type(field), pointer :: field_info
               !! include basis type, nubmer of local basis, etc.
      integer :: der_x_trial, der_y_trial, trial_basis_index, &
                 der_x_test, der_y_test, test_basis_index
      procedure(cofunc_2d), pointer, nopass, public :: cofunc => null()
      procedure(local_basis_func), pointer, nopass, public :: lbfunc => null()

   contains
      private
      procedure, public :: init_funcs
   end type func_to_quad_2d
contains

   subroutine init_funcs(self, cof_2d, lb_func, &
                         field_info, vertices, &
                         der_x_trial, der_y_trial, trial_basis_index, &
                         der_x_test, der_y_test, test_basis_index)
                  !! initialize the function to be integrated with input information.

      interface
         pure function local_b_func(x, y, local_vertices, basis_type, basis_index, der_x, der_y) result(func)
            import :: wp
            real(wp), intent(in) :: x, y
            real(wp), intent(in) :: local_vertices(2, 3)
            integer, intent(in) :: basis_type, basis_index, der_x, der_y
            real(wp) :: func
         end function local_b_func

         pure function cofunc_1d(x)
            import :: wp
            real(wp), intent(in) :: x
            real(wp) :: cofunc_1d
         end function cofunc_1d

         pure function cofunc_2d(x, y)
            import :: wp
            real(wp), intent(in) :: x, y
            real(wp) :: cofunc_2d
         end function cofunc_2d

         pure function cofunc_2d_t(x, y, t)
            import :: wp
            real(wp), intent(in) :: x, y, t
            real(wp) :: cofunc_2d_t
         end function cofunc_2d_t
      end interface

      class(func_to_quad_2d), intent(inout) :: self
      procedure(cofunc_2d), pointer, intent(in) :: cof_2d
      procedure(local_b_func), pointer, intent(in) ::lb_func
      class(field), target, intent(in) :: field_info
      integer, intent(in) :: der_x_trial, der_y_trial, trial_basis_index
      integer, intent(in) :: der_x_test, der_y_test, test_basis_index
      real(wp), intent(in) :: vertices(2, 3)

      self%cofunc => cof_2d
      self%lbfunc => lb_func
      self%field_info => field_info
      self%vertices = vertices
      self%der_x_trial = der_x_trial
      self%der_y_trial = der_y_trial
      self%trial_basis_index = trial_basis_index
      self%der_x_test = der_x_test
      self%der_y_test = der_y_test
      self%test_basis_index = test_basis_index

   end subroutine init_funcs

   real(wp) function afunc(self, x, y)
      !! Alogrithm I-3,  integrate function
      class(integration_class_2d), intent(inout) :: self
      real(wp), intent(in) :: x, y
      select type (self)
      class is (func_to_quad_2d)
         afunc = self%cofunc(x, y)* &
                 self%lbfunc(x, y, self%vertices, &
                             self%field_info%trial_basis_type, &
                             self%trial_basis_index, self%der_x_trial, &
                             self%der_y_trial)* &
                 self%lbfunc(x, y, self%vertices, &
                             self%field_info%test_basis_type, &
                             self%test_basis_index, self%der_x_test, &
                             self%der_y_test)
      class default
         error stop "Undefined Class"
      end select
   end function afunc

   real(wp) function bfunc(self, x, y)
      !! Alogrithm II-3, integrate function
      class(integration_class_2d), intent(inout) :: self
      real(wp), intent(in) :: x, y
      select type (self)
      class is (func_to_quad_2d)
         bfunc = self%cofunc(x, y)* &
                 self%lbfunc(x, y, self%vertices, &
                             self%field_info%test_basis_type, &
                             self%test_basis_index, self%der_x_test, &
                             self%der_y_test)
      class default
         error stop "Undefined Class"
      end select
   end function bfunc

end module mod_gauss_quad_2D
