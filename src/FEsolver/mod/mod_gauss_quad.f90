module mod_gauss_quad
   use quadrature_module, wp => quadrature_wp
   use mod_basis_func
   use mod_co_func
   implicit none

   private
   public :: quad_1D, quad_2D

   type, extends(integration_class_1d) :: quad_1D
      private
      class(func), pointer :: co_func
            !! The coefficient function of the integrand.
      class(local_basis), pointer :: l_basis
            !! local basis function
      real(wp) :: vertices(2)
            !! vertices: start[vertice(1)] and end[vertice(2)] of the local partition
      integer :: basis_type, basis_index, derivative_degree
            !! basis_type: the type of the FE.
            !! basis_type=101:1D linear FE.
            !! basis_index: the index of basis function to specify
            !!              which basis function we want to use.
            !! derivative_degree:the derivative degree of the FE basis function.
   contains
      procedure :: func_init => func_init_1D
      procedure :: quad_func => quad_func_1D
   end type quad_1D

   type, extends(integration_class_2d) :: quad_2D

   end type quad_2D

contains
   function quad_func_1D(self, x) result(f)
      class(quad_1D), intent(inout) :: self
      real(wp), intent(in) :: x
      real(wp) :: f

      f = self%l_basis%basis_f(x, &
                               self%vertices, &
                               self%basis_type, &
                               self%basis_index, &
                               self%derivative_degree)* &
          self%co_func%func_pass(x)

   end function quad_func_1D

   subroutine func_init_1D(self, co_func, l_basis, &
                           vertices, basis_type, &
                           basis_index, derivative_degree)
      class(quad_1D), intent(inout) :: self
      class(func), target :: co_func
      class(local_basis), target :: l_basis
      real(wp), intent(in) :: vertices(2)
      integer, intent(in) :: basis_type, basis_index, derivative_degree

      self%co_func => co_func
      self%l_basis => l_basis
      self%vertices = vertices
      self%basis_type = basis_type
      self%basis_index = basis_index
      self%derivative_degree = derivative_degree

   end subroutine func_init_1D
end module mod_gauss_quad
